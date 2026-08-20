#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
set.seed(20260809)

cohorts <- c("GSE75010_BIOBANK", "GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942")
out <- file.path(root, "results", "01_phase1b")
cohort_dir <- file.path(out, "cohort_DE")
meta_dir <- file.path(out, "meta")
robust_dir <- file.path(out, "robustness")
qc_dir <- file.path(out, "qc")
fig_dir <- file.path(out, "figures")
interim <- file.path(root, "data", "interim", "phase1b")
for (d in c(cohort_dir, meta_dir, robust_dir, qc_dir, fig_dir, interim)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "", quote = TRUE, fileEncoding = "UTF-8")
}

ga_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  weekly <- grepl("^[0-9]+w\\+[0-9]+d$", x)
  if (any(weekly)) {
    bits <- strsplit(x[weekly], "w\\+|d")
    out[weekly] <- vapply(bits, function(z) as.numeric(z[[1]]) + as.numeric(z[[2]]) / 7, numeric(1))
  }
  out
}

row_variance <- function(x) {
  means <- rowMeans(x)
  rowSums((x - means)^2) / (ncol(x) - 1)
}

max_vif <- function(design) {
  if (ncol(design) <= 2) return(1)
  values <- numeric(ncol(design) - 1)
  for (j in 2:ncol(design)) {
    y <- design[, j]
    others <- design[, -j, drop = FALSE]
    fit <- lm.fit(others, y)
    rss <- sum(fit$residuals^2)
    tss <- sum((y - mean(y))^2)
    r2 <- if (tss > 0) 1 - rss / tss else 1
    values[j - 1] <- 1 / max(1 - r2, .Machine$double.eps)
  }
  max(values)
}

read_gene_matrix <- function(path) {
  tab <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  genes <- as.character(tab[[1]])
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  stopifnot(!anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)))
  mat
}

freeze <- read.csv(file.path(root, "results", "01_phase1a", "bulk_sample_freeze.csv"), check.names = FALSE, stringsAsFactors = FALSE)
registry <- read.csv(file.path(root, "results", "01_phase1a1", "formal_phase1b_matrix_registry.csv"), check.names = FALSE, stringsAsFactors = FALSE)

matrices <- list()
metadata <- list()
sources <- setNames(registry$source, registry$dataset)
matrix_paths <- setNames(registry$gene_matrix_path, registry$dataset)
for (cohort in cohorts) {
  matrices[[cohort]] <- read_gene_matrix(file.path(root, matrix_paths[[cohort]]))
  md <- freeze[freeze$dataset == cohort & freeze$include_phase1b == "YES" & freeze$`PE/control` %in% c("PE", "CONTROL"), , drop = FALSE]
  md <- md[match(intersect(colnames(matrices[[cohort]]), md$`GSM/sample ID`), md$`GSM/sample ID`), , drop = FALSE]
  stopifnot(!anyNA(md$`GSM/sample ID`), all(md$`GSM/sample ID` %in% colnames(matrices[[cohort]])))
  metadata[[cohort]] <- md
}

make_model_data <- function(md) {
  data.frame(
    disease = factor(md$`PE/control`, levels = c("CONTROL", "PE")),
    GA_c = ga_numeric(md$GA) - mean(ga_numeric(md$GA), na.rm = TRUE),
    fetal_sex = factor(md$`fetal sex`, levels = c("F", "M")),
    batch = factor(md$batch, levels = c("A", "B")),
    labor = factor(md$labor, levels = c("spontaneous", "induced")),
    check.names = FALSE
  )
}

fit_limma_model <- function(cohort, formula_text, model_id) {
  md <- metadata[[cohort]]
  samples <- md$`GSM/sample ID`
  expr <- matrices[[cohort]][, samples, drop = FALSE]
  dat <- make_model_data(md)
  formula <- as.formula(formula_text)
  needed <- all.vars(formula)
  complete <- complete.cases(dat[, needed, drop = FALSE])
  dat <- dat[complete, , drop = FALSE]
  expr <- expr[, complete, drop = FALSE]
  design <- model.matrix(formula, dat)
  stopifnot(qr(design)$rank == ncol(design), "diseasePE" %in% colnames(design), nrow(design) > ncol(design))
  fit <- lmFit(expr, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  j <- match("diseasePE", colnames(design))
  effect <- fit$coefficients[, j]
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, j]
  p <- fit$p.value[, j]
  de <- data.frame(
    gene = rownames(expr), log2FC = effect, SE = se, t_statistic = fit$t[, j],
    raw_P = p, BH_FDR = p.adjust(p, method = "BH"),
    n_PE = sum(dat$disease == "PE"), n_control = sum(dat$disease == "CONTROL"),
    cohort = cohort, model_id = model_id, model_formula = formula_text,
    matrix_path = matrix_paths[[cohort]], mapping_version = "HGNC_2026-08-09_faaeb6ae1e2a596b",
    source = paste0(sources[[cohort]], "|results/01_phase1a/bulk_sample_freeze.csv"),
    stringsAsFactors = FALSE
  )
  categorical <- intersect(c("batch", "fetal_sex", "labor"), needed)
  zero_cells <- character()
  for (v in categorical) {
    xt <- table(dat$disease, dat[[v]])
    if (any(xt == 0)) zero_cells <- c(zero_cells, paste(v, paste(which(xt == 0, arr.ind = TRUE), collapse = ":"), sep = "="))
  }
  diagnostic <- data.frame(
    cohort = cohort, model_id = model_id, model_formula = formula_text,
    analysis_role = if (grepl("PRIMARY", model_id)) "PRIMARY_META_ESTIMAND" else "BIOLOGICAL_COVARIATE_SENSITIVITY",
    input_sample_n = length(samples), complete_case_n = nrow(design), n_PE = sum(dat$disease == "PE"), n_control = sum(dat$disease == "CONTROL"),
    design_column_n = ncol(design), design_rank = qr(design)$rank, residual_df = nrow(design) - ncol(design),
    condition_number = kappa(design), max_VIF = max_vif(design),
    disease_categorical_zero_cells = if (length(zero_cells)) paste(zero_cells, collapse = ";") else "NONE",
    coefficient = "diseasePE", limma_trend = TRUE, limma_robust = TRUE,
    estimability = "FULL_RANK_ESTIMABLE", sample_exclusions = if (all(complete)) "NONE" else paste(samples[!complete], collapse = ";"),
    source = paste0(sources[[cohort]], "|config/phase1b_analysis.json"), stringsAsFactors = FALSE
  )
  list(de = de, diagnostic = diagnostic, design = design)
}

primary_formulas <- c(
  GSE75010_BIOBANK = "~ disease", GSE30186 = "~ disease", GSE10588 = "~ disease",
  GSE24129 = "~ disease", GSE25906 = "~ disease + batch", GSE43942 = "~ disease"
)
primary <- list()
diagnostics <- list()
for (cohort in cohorts) {
  message("Fitting primary model: ", cohort)
  fit <- fit_limma_model(cohort, primary_formulas[[cohort]], paste0(cohort, "_PRIMARY"))
  primary[[cohort]] <- fit$de
  diagnostics[[length(diagnostics) + 1L]] <- fit$diagnostic
  write_csv(fit$de, file.path(cohort_dir, paste0(cohort, "_DE.csv")))
}

sensitivity_specs <- list(
  GSE75010_BIOBANK_GA = c("GSE75010_BIOBANK", "~ disease + GA_c"),
  GSE75010_BIOBANK_GA_SEX = c("GSE75010_BIOBANK", "~ disease + GA_c + fetal_sex"),
  GSE25906_GA = c("GSE25906", "~ disease + batch + GA_c"),
  GSE25906_GA_SEX = c("GSE25906", "~ disease + batch + GA_c + fetal_sex"),
  GSE25906_GA_SEX_LABOR = c("GSE25906", "~ disease + batch + GA_c + fetal_sex + labor")
)
sens <- list()
for (model_id in names(sensitivity_specs)) {
  spec <- sensitivity_specs[[model_id]]
  message("Fitting sensitivity model: ", model_id)
  fit <- fit_limma_model(spec[[1]], spec[[2]], model_id)
  if (fit$diagnostic$max_VIF >= 5 || fit$diagnostic$disease_categorical_zero_cells != "NONE") stop("Sensitivity design failed stability lock: ", model_id)
  sens[[model_id]] <- fit$de
  diagnostics[[length(diagnostics) + 1L]] <- fit$diagnostic
  write_csv(fit$de, file.path(interim, paste0(model_id, "_DE.csv")))
}

# Gene availability is explicit and long-form; an absent platform entry is never a biological zero.
gene_union <- sort(unique(unlist(lapply(matrices, rownames), use.names = FALSE)))
availability <- do.call(rbind, lapply(cohorts, function(cohort) {
  present <- gene_union %in% rownames(matrices[[cohort]])
  data.frame(
    gene = gene_union, cohort = cohort,
    platform_available = ifelse(present, "YES", "NO"),
    mapping_status = ifelse(present, "MAPPED_UNAMBIGUOUS", "NOT_AVAILABLE_OR_NOT_UNAMBIGUOUS"),
    estimable_primary = ifelse(present, "YES", "NO"),
    absence_interpretation = ifelse(present, "NOT_APPLICABLE", "PLATFORM_OR_MAPPING_ABSENCE_NOT_BIOLOGICAL_NEGATIVE"),
    matrix_path = matrix_paths[[cohort]], mapping_version = "HGNC_2026-08-09_faaeb6ae1e2a596b",
    source = paste0(sources[[cohort]], "|results/01_phase1a/gene_mapping_registry.csv"),
    stringsAsFactors = FALSE
  )
}))
write_csv(availability, file.path(qc_dir, "phase1b_gene_availability.csv"))

# Vectorized intercept-only REML with conservative modified Knapp-Hartung t inference.
reml_score <- function(y, v, tau2) {
  w <- 1 / (v + tau2)
  sw <- rowSums(w)
  mu <- rowSums(w * y) / sw
  residual <- y - mu
  0.5 * (rowSums((w^2) * (residual^2)) - (sw - rowSums(w^2) / sw))
}

meta_block <- function(y, se) {
  v <- se^2
  n <- nrow(y); k <- ncol(y)
  score0 <- reml_score(y, v, rep(0, n))
  active <- is.finite(score0) & score0 > 0
  tau2 <- rep(0, n)
  if (any(active)) {
    lo <- rep(0, n)
    empirical <- apply(y, 1, var)
    hi <- pmax(empirical, rowMeans(v), 1e-4) * 10 + 1e-4
    for (iter in 1:30) {
      shi <- reml_score(y, v, hi)
      need <- active & is.finite(shi) & shi > 0
      if (!any(need)) break
      hi[need] <- hi[need] * 4
    }
    for (iter in 1:60) {
      mid <- (lo + hi) / 2
      smid <- reml_score(y, v, mid)
      move_lo <- active & is.finite(smid) & smid > 0
      lo[move_lo] <- mid[move_lo]
      hi[active & !move_lo] <- mid[active & !move_lo]
    }
    tau2[active] <- (lo[active] + hi[active]) / 2
  }
  w <- 1 / (v + tau2)
  sw <- rowSums(w)
  mu <- rowSums(w * y) / sw
  qe_random <- rowSums(w * (y - mu)^2)
  kh_scale <- pmax(1, qe_random / (k - 1))
  pooled_se <- sqrt(kh_scale / sw)
  stat <- mu / pooled_se
  df <- k - 1
  p <- 2 * pt(-abs(stat), df = df)
  crit <- qt(0.975, df = df)
  w0 <- 1 / v
  mu0 <- rowSums(w0 * y) / rowSums(w0)
  q <- rowSums(w0 * (y - mu0)^2)
  i2 <- ifelse(q > 0, pmax(0, (q - (k - 1)) / q) * 100, 0)
  data.frame(
    pooled_effect = mu, pooled_SE = pooled_se, CI_lower = mu - crit * pooled_se, CI_upper = mu + crit * pooled_se,
    test_statistic = stat, test_df = df, raw_meta_P = p, tau2 = tau2, I2 = i2,
    Cochran_Q = q, Q_P = pchisq(q, df = k - 1, lower.tail = FALSE), stringsAsFactors = FALSE
  )
}

meta_matrix <- function(y, se, minimum_k) {
  ok <- is.finite(y) & is.finite(se) & se > 0
  codes <- apply(ok, 1, function(z) paste(as.integer(z), collapse = ""))
  result <- data.frame(
    pooled_effect = rep(NA_real_, nrow(y)), pooled_SE = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
    test_statistic = NA_real_, test_df = NA_real_, raw_meta_P = NA_real_, tau2 = NA_real_, I2 = NA_real_,
    Cochran_Q = NA_real_, Q_P = NA_real_, stringsAsFactors = FALSE
  )
  for (code in unique(codes)) {
    cols <- strsplit(code, "", fixed = TRUE)[[1]] == "1"
    if (sum(cols) < minimum_k) next
    rows <- which(codes == code)
    block <- meta_block(y[rows, cols, drop = FALSE], se[rows, cols, drop = FALSE])
    result[rows, ] <- block
  }
  result$k_cohorts <- rowSums(ok)
  result
}

effect_matrix <- matrix(NA_real_, nrow = length(gene_union), ncol = length(cohorts), dimnames = list(gene_union, cohorts))
se_matrix <- effect_matrix
for (cohort in cohorts) {
  idx <- match(primary[[cohort]]$gene, gene_union)
  effect_matrix[idx, cohort] <- primary[[cohort]]$log2FC
  se_matrix[idx, cohort] <- primary[[cohort]]$SE
}
eligible <- rowSums(is.finite(effect_matrix) & is.finite(se_matrix) & se_matrix > 0) >= 4
eligible_genes <- gene_union[eligible]
y <- effect_matrix[eligible, , drop = FALSE]
se <- se_matrix[eligible, , drop = FALSE]
message("Primary meta-analysis eligible genes: ", nrow(y))
meta_fit <- meta_matrix(y, se, minimum_k = 4)

# Validate the custom vectorized REML tau2/effect against metafor on alphabetically selected genes.
validation_idx <- seq_len(min(40, nrow(y)))
tau_diff <- effect_diff <- numeric(length(validation_idx))
for (ii in seq_along(validation_idx)) {
  i <- validation_idx[[ii]]
  keep <- is.finite(y[i, ]) & is.finite(se[i, ]) & se[i, ] > 0
  mf <- suppressWarnings(rma.uni(yi = y[i, keep], sei = se[i, keep], method = "REML"))
  tau_diff[ii] <- abs(as.numeric(mf$tau2) - meta_fit$tau2[i])
  effect_diff[ii] <- abs(as.numeric(mf$beta[1]) - meta_fit$pooled_effect[i])
}
message("REML validation max tau2 diff: ", format(max(tau_diff), scientific = TRUE),
        "; max pooled-effect diff: ", format(max(effect_diff), scientific = TRUE))
stopifnot(max(tau_diff) < 1e-4, max(effect_diff) < 1e-4)

same_direction <- numeric(nrow(y)); n_pos <- rowSums(y > 0, na.rm = TRUE); n_neg <- rowSums(y < 0, na.rm = TRUE)
for (i in seq_len(nrow(y))) {
  same_direction[i] <- if (meta_fit$pooled_effect[i] >= 0) n_pos[i] / meta_fit$k_cohorts[i] else n_neg[i] / meta_fit$k_cohorts[i]
}
meta_fit$BH_FDR <- p.adjust(meta_fit$raw_meta_P, method = "BH")

# LOCO p-values are BH-adjusted over the full valid gene family in each omission.
loco_all <- list()
for (omit in cohorts) {
  message("LOCO meta-analysis omitting: ", omit)
  yy <- y; ss <- se
  yy[, omit] <- NA_real_; ss[, omit] <- NA_real_
  lf <- meta_matrix(yy, ss, minimum_k = 3)
  lf$BH_FDR <- NA_real_
  valid <- is.finite(lf$raw_meta_P)
  lf$BH_FDR[valid] <- p.adjust(lf$raw_meta_P[valid], method = "BH")
  loco_all[[omit]] <- lf
}

loco_same_all <- logical(nrow(y)); loco_fdr_prop <- numeric(nrow(y)); loco_valid_n <- integer(nrow(y))
for (i in seq_len(nrow(y))) {
  effects <- vapply(loco_all, function(z) z$pooled_effect[i], numeric(1))
  fdrs <- vapply(loco_all, function(z) z$BH_FDR[i], numeric(1))
  valid <- is.finite(effects) & is.finite(fdrs)
  loco_valid_n[i] <- sum(valid)
  loco_same_all[i] <- all(sign(effects[valid]) == sign(meta_fit$pooled_effect[i]))
  loco_fdr_prop[i] <- if (any(valid)) mean(fdrs[valid] < 0.10) else 0
}

prelim_robust <- meta_fit$BH_FDR < 0.05 & same_direction >= 0.75
stable <- prelim_robust & meta_fit$I2 <= 60 & loco_same_all & loco_fdr_prop >= 0.80
robust_heterogeneous <- prelim_robust & meta_fit$I2 > 60 & loco_same_all & loco_fdr_prop >= 0.80
category <- rep("UNSTABLE", nrow(y))
category[meta_fit$BH_FDR < 0.05 & !(stable | robust_heterogeneous)] <- "COHORT_SPECIFIC"
category[meta_fit$BH_FDR >= 0.05 & same_direction >= 0.75] <- "DIRECTION_CONSISTENT_NON_SIGNIFICANT"
category[robust_heterogeneous] <- "ROBUST_BUT_HETEROGENEOUS"
category[stable] <- "STABLE"

meta_table <- data.frame(
  gene = eligible_genes, k_cohorts = meta_fit$k_cohorts,
  pooled_log2FC = meta_fit$pooled_effect, pooled_SE = meta_fit$pooled_SE,
  CI_lower = meta_fit$CI_lower, CI_upper = meta_fit$CI_upper,
  test_statistic = meta_fit$test_statistic, test_df = meta_fit$test_df,
  raw_meta_P = meta_fit$raw_meta_P, BH_FDR = meta_fit$BH_FDR,
  tau2 = meta_fit$tau2, I2 = meta_fit$I2, Cochran_Q = meta_fit$Cochran_Q, Q_P = meta_fit$Q_P,
  direction_consistency = same_direction, n_positive = n_pos, n_negative = n_neg,
  valid_LOCO_n = loco_valid_n, all_LOCO_same_direction = loco_same_all,
  LOCO_FDR_lt_0_10_proportion = loco_fdr_prop, category = category,
  meta_method = "REML_MODIFIED_KNHA_T", mapping_version = "HGNC_2026-08-09_faaeb6ae1e2a596b",
  source = "six frozen cohort DE tables|config/phase1b_analysis.json",
  stringsAsFactors = FALSE
)
for (cohort in cohorts) {
  meta_table[[paste0(cohort, "_log2FC")]] <- y[, cohort]
  meta_table[[paste0(cohort, "_SE")]] <- se[, cohort]
}
meta_table$cohort_effects <- apply(y, 1, function(z) paste(paste(cohorts[is.finite(z)], format(z[is.finite(z)], digits = 6), sep = "="), collapse = ";"))
meta_table <- meta_table[order(meta_table$BH_FDR, -abs(meta_table$pooled_log2FC), meta_table$gene), ]
write_csv(meta_table, file.path(meta_dir, "pe_gene_meta_analysis.csv"))
write_csv(meta_table[meta_table$category == "STABLE", ], file.path(meta_dir, "stable_pe_genes.csv"))
write_csv(meta_table[meta_table$I2 > 60, ], file.path(meta_dir, "heterogeneous_pe_genes.csv"))

direction_table <- meta_table[, c("gene", "k_cohorts", "pooled_log2FC", "BH_FDR", "direction_consistency", "n_positive", "n_negative", "category", paste0(cohorts, "_log2FC"), "source")]
write_csv(direction_table, file.path(meta_dir, "direction_consistency.csv"))

# Persist robust-candidate LOCO rows only; BH values came from the full LOCO families above.
candidate_genes <- eligible_genes[prelim_robust]
loco_rows <- list()
for (omit in cohorts) {
  lf <- loco_all[[omit]]
  idx <- which(prelim_robust & is.finite(lf$pooled_effect))
  if (!length(idx)) next
  loco_rows[[length(loco_rows) + 1L]] <- data.frame(
    gene = eligible_genes[idx], omitted_cohort = omit, remaining_cohort_n = lf$k_cohorts[idx],
    pooled_log2FC = lf$pooled_effect[idx], pooled_SE = lf$pooled_SE[idx], raw_meta_P = lf$raw_meta_P[idx],
    BH_FDR = lf$BH_FDR[idx], I2 = lf$I2[idx], direction = ifelse(lf$pooled_effect[idx] >= 0, "UP", "DOWN"),
    full_pooled_log2FC = meta_fit$pooled_effect[idx], same_direction_as_full = sign(lf$pooled_effect[idx]) == sign(meta_fit$pooled_effect[idx]),
    full_category = category[idx], BH_family = "all genes valid with >=3 remaining cohorts",
    source = "six frozen cohort DE tables|config/phase1b_analysis.json", stringsAsFactors = FALSE
  )
}
loco_table <- if (length(loco_rows)) do.call(rbind, loco_rows) else data.frame(
  gene = character(), omitted_cohort = character(), remaining_cohort_n = integer(), pooled_log2FC = numeric(), pooled_SE = numeric(),
  raw_meta_P = numeric(), BH_FDR = numeric(), I2 = numeric(), direction = character(), full_pooled_log2FC = numeric(),
  same_direction_as_full = logical(), full_category = character(), BH_family = character(), source = character()
)
write_csv(loco_table, file.path(robust_dir, "leave_one_cohort_out.csv"))

# Unadjusted Hedges' g by cohort: platform-dynamic-range sensitivity only.
g_matrix <- matrix(NA_real_, nrow = length(gene_union), ncol = length(cohorts), dimnames = list(gene_union, cohorts))
g_se_matrix <- g_matrix
for (cohort in cohorts) {
  md <- metadata[[cohort]]
  expr <- matrices[[cohort]][, md$`GSM/sample ID`, drop = FALSE]
  pe <- md$`PE/control` == "PE"; control <- md$`PE/control` == "CONTROL"
  n1 <- sum(pe); n0 <- sum(control); total <- n1 + n0
  pooled_sd <- sqrt(((n1 - 1) * row_variance(expr[, pe, drop = FALSE]) + (n0 - 1) * row_variance(expr[, control, drop = FALSE])) / (total - 2))
  d <- (rowMeans(expr[, pe, drop = FALSE]) - rowMeans(expr[, control, drop = FALSE])) / pooled_sd
  correction <- 1 - 3 / (4 * total - 9)
  g <- correction * d
  vg <- total / (n1 * n0) + g^2 / (2 * (total - 2))
  invalid <- !is.finite(g) | !is.finite(vg) | vg <= 0
  g[invalid] <- NA_real_; vg[invalid] <- NA_real_
  idx <- match(rownames(expr), gene_union)
  g_matrix[idx, cohort] <- g
  g_se_matrix[idx, cohort] <- sqrt(vg)
}
gy <- g_matrix[eligible, , drop = FALSE]; gse <- g_se_matrix[eligible, , drop = FALSE]
g_meta <- meta_matrix(gy, gse, minimum_k = 4)
g_meta$BH_FDR <- NA_real_
gvalid <- is.finite(g_meta$raw_meta_P)
g_meta$BH_FDR[gvalid] <- p.adjust(g_meta$raw_meta_P[gvalid], method = "BH")
g_pos <- rowSums(gy > 0, na.rm = TRUE); g_neg <- rowSums(gy < 0, na.rm = TRUE)
g_direction_consistency <- ifelse(g_meta$pooled_effect >= 0, g_pos / g_meta$k_cohorts, g_neg / g_meta$k_cohorts)
std_table <- data.frame(
  gene = eligible_genes, standardized_k_cohorts = g_meta$k_cohorts,
  pooled_Hedges_g = g_meta$pooled_effect, pooled_SE = g_meta$pooled_SE,
  CI_lower = g_meta$CI_lower, CI_upper = g_meta$CI_upper, test_statistic = g_meta$test_statistic,
  test_df = g_meta$test_df, raw_meta_P = g_meta$raw_meta_P, BH_FDR = g_meta$BH_FDR,
  tau2 = g_meta$tau2, I2 = g_meta$I2, Cochran_Q = g_meta$Cochran_Q, Q_P = g_meta$Q_P,
  standardized_direction_consistency = g_direction_consistency, standardized_n_positive = g_pos, standardized_n_negative = g_neg,
  primary_pooled_log2FC = meta_fit$pooled_effect, primary_BH_FDR = meta_fit$BH_FDR,
  primary_category = category, direction_matches_primary = sign(g_meta$pooled_effect) == sign(meta_fit$pooled_effect),
  is_primary_STABLE = stable, standardized_significant_same_direction_as_primary = g_meta$BH_FDR < 0.05 & sign(g_meta$pooled_effect) == sign(meta_fit$pooled_effect),
  interpretation = "PLATFORM_DYNAMIC_RANGE_SENSITIVITY_NOT_COVARIATE_ADJUSTED",
  source = "frozen formal matrices|config/phase1b_analysis.json", stringsAsFactors = FALSE
)
std_table$primary_rank <- rank(std_table$primary_BH_FDR, ties.method = "min", na.last = "keep")
std_table$standardized_rank <- rank(std_table$BH_FDR, ties.method = "min", na.last = "keep")
write_csv(std_table[order(std_table$BH_FDR, std_table$gene), ], file.path(robust_dir, "standardized_effect_sensitivity.csv"))

# Covariate sensitivity is restricted to primary stable genes and never alters membership.
stable_genes <- eligible_genes[stable]
de_lookup <- function(tab, genes) tab[match(genes, tab$gene), c("gene", "log2FC", "SE"), drop = FALSE]
cov_rows <- list(); flag_sets <- list()
add_cov_comparison <- function(cohort, comparison, reference_id, adjusted_id, flag_type) {
  genes <- intersect(stable_genes, primary[[cohort]]$gene)
  if (!length(genes)) return(NULL)
  ptab <- de_lookup(primary[[cohort]], genes)
  rtab <- if (reference_id == "PRIMARY") ptab else de_lookup(sens[[reference_id]], genes)
  atab <- de_lookup(sens[[adjusted_id]], genes)
  attenuation <- abs(atab$log2FC) / abs(rtab$log2FC)
  incremental <- abs(atab$log2FC - rtab$log2FC) / pmax(abs(rtab$log2FC), 0.05)
  flip <- sign(atab$log2FC) != sign(rtab$log2FC)
  sensitive <- if (flag_type == "GA") flip | attenuation < 0.5 else flip | incremental > 0.25
  flag_sets[[paste(cohort, flag_type, sep = "_")]] <<- setNames(sensitive, genes)
  correlation <- if (length(genes) >= 3) cor(ptab$log2FC, atab$log2FC, use = "complete.obs") else NA_real_
  cov_rows[[length(cov_rows) + 1L]] <<- data.frame(
    gene = genes, cohort = cohort, comparison = comparison, reference_model = reference_id, adjusted_model = adjusted_id,
    primary_log2FC = ptab$log2FC, reference_log2FC = rtab$log2FC, adjusted_log2FC = atab$log2FC,
    sign_concordant_with_primary = sign(atab$log2FC) == sign(ptab$log2FC),
    attenuation_ratio_vs_primary = abs(atab$log2FC) / abs(ptab$log2FC), attenuation_ratio_vs_reference = attenuation,
    incremental_relative_change = incremental, direction_flip_vs_reference = flip,
    GA_SENSITIVE = flag_type == "GA" & sensitive, SEX_SENSITIVE = flag_type == "SEX" & sensitive,
    LABOR_SENSITIVE = flag_type == "LABOR" & sensitive,
    effect_correlation_with_primary_across_stable_genes = correlation, stable_gene_n_in_cohort = length(genes),
    interpretation = "ANNOTATION_ONLY_NO_GENE_REMOVAL", source = paste0(sources[[cohort]], "|config/phase1b_analysis.json"),
    stringsAsFactors = FALSE
  )
}
invisible(add_cov_comparison("GSE75010_BIOBANK", "GA_COMPONENT", "PRIMARY", "GSE75010_BIOBANK_GA", "GA"))
invisible(add_cov_comparison("GSE75010_BIOBANK", "SEX_COMPONENT", "GSE75010_BIOBANK_GA", "GSE75010_BIOBANK_GA_SEX", "SEX"))
invisible(add_cov_comparison("GSE25906", "GA_COMPONENT", "PRIMARY", "GSE25906_GA", "GA"))
invisible(add_cov_comparison("GSE25906", "SEX_COMPONENT", "GSE25906_GA", "GSE25906_GA_SEX", "SEX"))
invisible(add_cov_comparison("GSE25906", "LABOR_COMPONENT", "GSE25906_GA_SEX", "GSE25906_GA_SEX_LABOR", "LABOR"))
cov_table <- if (length(cov_rows)) do.call(rbind, cov_rows) else data.frame(
  gene = character(), cohort = character(), comparison = character(), reference_model = character(), adjusted_model = character(),
  primary_log2FC = numeric(), reference_log2FC = numeric(), adjusted_log2FC = numeric(), sign_concordant_with_primary = logical(),
  attenuation_ratio_vs_primary = numeric(), attenuation_ratio_vs_reference = numeric(), incremental_relative_change = numeric(), direction_flip_vs_reference = logical(),
  GA_SENSITIVE = logical(), SEX_SENSITIVE = logical(), LABOR_SENSITIVE = logical(),
  effect_correlation_with_primary_across_stable_genes = numeric(), stable_gene_n_in_cohort = integer(),
  interpretation = character(), source = character()
)
write_csv(cov_table, file.path(robust_dir, "covariate_sensitivity.csv"))

# Model diagnostics and meta-engine validation.
diag_table <- do.call(rbind, diagnostics)
meta_diag <- data.frame(
  cohort = "CROSS_COHORT", model_id = "REML_ENGINE_VALIDATION", model_formula = "intercept-only random effects",
  analysis_role = "TECHNICAL_VALIDATION", input_sample_n = nrow(y), complete_case_n = length(validation_idx),
  n_PE = NA, n_control = NA, design_column_n = NA, design_rank = NA, residual_df = NA,
  condition_number = NA, max_VIF = NA, disease_categorical_zero_cells = "NOT_APPLICABLE",
  coefficient = "pooled_effect", limma_trend = NA, limma_robust = NA,
  estimability = paste0("MAX_TAU2_DIFF_METAFOR=", format(max(tau_diff), scientific = TRUE), ";MAX_EFFECT_DIFF=", format(max(effect_diff), scientific = TRUE)),
  sample_exclusions = "NONE", source = "metafor 5.0.1 validation of first 40 alphabetical eligible genes",
  stringsAsFactors = FALSE
)
diag_table <- rbind(diag_table, meta_diag)
write_csv(diag_table, file.path(qc_dir, "phase1b_model_diagnostics.csv"))

# Result-aware risk flags do not alter thresholds or membership.
stable_n <- sum(stable); hetero_robust_n <- sum(robust_heterogeneous)
risk_rows <- data.frame(
  risk_id = c("P1B-R01", "P1B-R02", "P1B-R03", "P1B-R04", "P1B-R05", "P1B-R06"),
  severity = c("HIGH", "HIGH", "MODERATE", "MODERATE", if (stable_n == 0) "CRITICAL" else "MODERATE", "MODERATE"),
  scope = c("GSE30186;GSE10588;GSE43942", "GSE25906", "GSE75010_BIOBANK", "CROSS_PLATFORM", "STABLE_SIGNATURE", "HETEROGENEITY"),
  risk = c(
    "Sample-level GA is unavailable in three core cohorts; the primary estimand is minimally adjusted and cannot remove that confounding.",
    "Batch is technically structured and mandatory in the GSE25906 primary model; biological covariates are sensitivity-only.",
    "Sample-level processing batch is unavailable for the BioBank cohort.",
    "Historical platforms have different dynamic ranges; log2 effects remain scale-dependent.",
    paste0("Stable-gene count under the frozen rule is ", stable_n, "."),
    paste0(hetero_robust_n, " genes meet all robustness criteria except I2 and remain non-STABLE.")
  ),
  mitigation = c(
    "Interpret cohort/meta effects as PE-associated rather than causally GA-independent; use available adjusted sensitivities.",
    "Retain all samples, include batch in primary and sensitivity models, and report labor/GA/sex sensitivities.",
    "Do not infer or impute batch; use LOCO and influence reporting.",
    "Use Hedges-g sensitivity and require direction/LOCO consistency; do not use a universal meta-log2FC cutoff.",
    "Do not relax thresholds; apply the preregistered Phase 1C gate.",
    "Label exploratory subtype/platform hypotheses; do not promote heterogeneous genes to STABLE."
  ),
  status = c("OPEN_INTERPRETIVE_RESTRICTION", "MODELED_RESTRICTION", "OPEN_INTERPRETIVE_RESTRICTION", "SENSITIVITY_COMPLETED", "FROZEN_RESULT", "FROZEN_RESULT"),
  source = c(
    "results/01_phase1a/bulk_sample_freeze.csv", "results/01_phase1b/qc/phase1b_model_diagnostics.csv",
    sources[["GSE75010_BIOBANK"]], "results/01_phase1b/robustness/standardized_effect_sensitivity.csv",
    "results/01_phase1b/meta/stable_pe_genes.csv", "results/01_phase1b/meta/heterogeneous_pe_genes.csv"
  ), stringsAsFactors = FALSE
)
write_csv(risk_rows, file.path(qc_dir, "phase1b_risk_flags.csv"))

# Analytical figures with algorithmic selection only.
stable_table <- meta_table[meta_table$category == "STABLE", ]
representative <- head(stable_table[order(stable_table$BH_FDR, stable_table$I2, stable_table$gene), ], 10)
plot_empty <- function(path, title, text) {
  png(path, width = 1400, height = 900, res = 150)
  plot.new(); title(main = title); text(0.5, 0.5, text, cex = 1.2); invisible(dev.off())
}

png(file.path(fig_dir, "A_meta_volcano.png"), width = 1600, height = 1200, res = 160)
cols <- ifelse(meta_table$category == "STABLE", "#B2182B", ifelse(meta_table$category == "ROBUST_BUT_HETEROGENEOUS", "#EF8A62", "#BDBDBD"))
plot(meta_table$pooled_log2FC, -log10(pmax(meta_table$BH_FDR, .Machine$double.xmin)), pch = 16, cex = 0.45, col = cols,
     xlab = "Random-effects pooled log2FC", ylab = "-log10(meta BH FDR)", main = "PE multi-cohort meta-analysis (algorithmic categories)")
abline(h = -log10(0.05), lty = 2, col = "#555555")
if (nrow(representative)) text(representative$pooled_log2FC, -log10(representative$BH_FDR), labels = representative$gene, pos = 3, cex = 0.65)
legend("topright", legend = c("STABLE", "ROBUST_BUT_HETEROGENEOUS", "Other"), col = c("#B2182B", "#EF8A62", "#BDBDBD"), pch = 16, bty = "n")
invisible(dev.off())

png(file.path(fig_dir, "F_heterogeneity_distribution.png"), width = 1400, height = 900, res = 150)
hist(meta_table$I2, breaks = seq(0, 100, by = 5), col = "#6BAED6", border = "white", xlab = "I2 (%)", main = "Cross-cohort heterogeneity distribution")
abline(v = 60, col = "#B2182B", lwd = 2, lty = 2); invisible(dev.off())

png(file.path(fig_dir, "E_standardized_effect_concordance.png"), width = 1300, height = 1100, res = 150)
valid_std <- is.finite(std_table$pooled_Hedges_g) & is.finite(std_table$primary_pooled_log2FC)
plot(std_table$primary_pooled_log2FC[valid_std], std_table$pooled_Hedges_g[valid_std], pch = 16, cex = 0.45,
     col = ifelse(std_table$is_primary_STABLE[valid_std], "#B2182B", "#BDBDBD"),
     xlab = "Primary pooled log2FC", ylab = "Pooled unadjusted Hedges' g", main = "Platform dynamic-range sensitivity")
abline(h = 0, v = 0, col = "#777777", lty = 3); invisible(dev.off())

heat_genes <- head(stable_table[order(stable_table$BH_FDR, stable_table$I2, stable_table$gene), "gene"], 50)
if (length(heat_genes)) {
  hm <- y[match(heat_genes, eligible_genes), , drop = FALSE]
  lim <- max(abs(hm), na.rm = TRUE)
  png(file.path(fig_dir, "B_cohort_effect_heatmap.png"), width = 1400, height = max(900, 24 * nrow(hm)), res = 150)
  par(mar = c(9, 9, 4, 2))
  image(seq_len(ncol(hm)), seq_len(nrow(hm)), t(hm[nrow(hm):1, , drop = FALSE]), col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "", main = "Cohort effects: top STABLE genes by meta FDR")
  axis(1, at = seq_len(ncol(hm)), labels = colnames(hm), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nrow(hm)), labels = rev(rownames(hm)), las = 2, cex.axis = 0.55)
  invisible(dev.off())
} else plot_empty(file.path(fig_dir, "B_cohort_effect_heatmap.png"), "Cohort effect heatmap", "No STABLE genes under the frozen rule")

forest_genes <- head(stable_table[order(stable_table$BH_FDR, stable_table$I2, stable_table$gene), "gene"], 6)
if (length(forest_genes)) {
  png(file.path(fig_dir, "A_representative_stable_forest.png"), width = 1800, height = 1200, res = 150)
  par(mfrow = c(2, 3), mar = c(4, 7, 3, 1))
  for (gene in forest_genes) {
    i <- match(gene, eligible_genes); keep <- is.finite(y[i, ])
    eff <- c(y[i, keep], meta_fit$pooled_effect[i]); serr <- c(se[i, keep], meta_fit$pooled_SE[i])
    labels <- c(cohorts[keep], "REML meta")
    yy <- rev(seq_along(eff)); limits <- range(eff - 2 * serr, eff + 2 * serr)
    plot(eff, yy, xlim = limits, ylim = c(0.5, length(eff) + 0.5), yaxt = "n", pch = c(rep(16, length(eff) - 1), 18), xlab = "log2FC", ylab = "", main = gene)
    segments(eff - 1.96 * serr, yy, eff + 1.96 * serr, yy); abline(v = 0, lty = 3)
    axis(2, at = yy, labels = labels, las = 2, cex.axis = 0.65)
  }
  invisible(dev.off())
} else plot_empty(file.path(fig_dir, "A_representative_stable_forest.png"), "Representative stable-gene forest plot", "No STABLE genes under the frozen rule")

loo_genes <- head(stable_table[order(stable_table$BH_FDR, stable_table$I2, stable_table$gene), "gene"], 15)
if (length(loo_genes)) {
  png(file.path(fig_dir, "D_leave_one_cohort_out_robustness.png"), width = 1600, height = 1000, res = 150)
  xlabels <- c("Full", paste0("omit ", cohorts))
  mat <- matrix(NA_real_, nrow = length(loo_genes), ncol = length(xlabels), dimnames = list(loo_genes, xlabels))
  for (i in seq_along(loo_genes)) {
    idx <- match(loo_genes[i], eligible_genes)
    mat[i, 1] <- meta_fit$pooled_effect[idx]
    mat[i, -1] <- vapply(loco_all, function(z) z$pooled_effect[idx], numeric(1))
  }
  matplot(t(mat), type = "l", lty = 1, lwd = 1.2, col = rainbow(nrow(mat)), xaxt = "n", xlab = "", ylab = "Pooled log2FC", main = "LOCO effects: top STABLE genes by meta FDR")
  axis(1, at = seq_along(xlabels), labels = xlabels, las = 2, cex.axis = 0.7); abline(h = 0, lty = 3)
  legend("topright", legend = rownames(mat), col = rainbow(nrow(mat)), lty = 1, cex = 0.55, bty = "n")
  invisible(dev.off())
} else plot_empty(file.path(fig_dir, "D_leave_one_cohort_out_robustness.png"), "LOCO robustness", "No STABLE genes under the frozen rule")

# Compact machine-readable summary for report construction.
std_support <- if (stable_n) mean(std_table$direction_matches_primary[match(stable_genes, std_table$gene)], na.rm = TRUE) else NA_real_
bio_omit_idx <- if (nrow(loco_table)) loco_table$omitted_cohort == "GSE75010_BIOBANK" & loco_table$full_category == "STABLE" else logical()
biobank_support <- if (any(bio_omit_idx)) mean(loco_table$same_direction_as_full[bio_omit_idx] & loco_table$BH_FDR[bio_omit_idx] < 0.10) else NA_real_
bio_all <- loco_all[["GSE75010_BIOBANK"]]
bio_valid <- is.finite(bio_all$pooled_effect) & is.finite(meta_fit$pooled_effect)
std_valid_summary <- is.finite(std_table$pooled_Hedges_g) & is.finite(std_table$primary_pooled_log2FC)
compare_full_sensitivity <- function(cohort, adjusted_id) {
  p <- primary[[cohort]]; a <- sens[[adjusted_id]]
  idx <- match(p$gene, a$gene); valid <- !is.na(idx)
  pe <- p$log2FC[valid]; ae <- a$log2FC[idx[valid]]
  c(correlation = cor(pe, ae), sign_agreement = mean(sign(pe) == sign(ae)),
    median_attenuation = median(abs(ae) / pmax(abs(pe), 1e-8)), direction_flip_n = sum(sign(pe) != sign(ae)))
}
bio_sens_summary <- compare_full_sensitivity("GSE75010_BIOBANK", "GSE75010_BIOBANK_GA_SEX")
g259_sens_summary <- compare_full_sensitivity("GSE25906", "GSE25906_GA_SEX")
g259_labor_summary <- compare_full_sensitivity("GSE25906", "GSE25906_GA_SEX_LABOR")
summary <- data.frame(
  metric = c(
    "gene_union_n", "estimable_ge4_n", "meta_fdr_lt_0_05_n", "stable_n", "robust_but_heterogeneous_n",
    "direction_consistent_non_significant_n", "cohort_specific_n", "unstable_n", "i2_gt_60_n", "minimum_meta_BH_FDR",
    "standardized_stable_direction_support", "standardized_all_direction_support", "standardized_effect_correlation", "standardized_rank_spearman",
    "biobank_omission_stable_support", "biobank_omission_all_sign_agreement", "biobank_omission_all_effect_correlation", "biobank_omission_meta_fdr_lt_0_05_n",
    "GSE75010_GA_sex_effect_correlation", "GSE75010_GA_sex_sign_agreement", "GSE75010_GA_sex_median_attenuation", "GSE75010_GA_sex_direction_flip_n",
    "GSE25906_GA_sex_effect_correlation", "GSE25906_GA_sex_sign_agreement", "GSE25906_GA_sex_median_attenuation", "GSE25906_GA_sex_direction_flip_n",
    "GSE25906_labor_effect_correlation", "GSE25906_labor_sign_agreement", "GSE25906_labor_median_attenuation", "GSE25906_labor_direction_flip_n",
    "meta_engine_max_tau2_diff", "meta_engine_max_effect_diff"
  ),
  value = c(
    length(gene_union), nrow(y), sum(meta_fit$BH_FDR < 0.05), stable_n, hetero_robust_n,
    sum(category == "DIRECTION_CONSISTENT_NON_SIGNIFICANT"), sum(category == "COHORT_SPECIFIC"), sum(category == "UNSTABLE"), sum(meta_fit$I2 > 60), min(meta_fit$BH_FDR),
    std_support, mean(std_table$direction_matches_primary[std_valid_summary]), cor(std_table$primary_pooled_log2FC[std_valid_summary], std_table$pooled_Hedges_g[std_valid_summary]), cor(std_table$primary_rank[std_valid_summary], std_table$standardized_rank[std_valid_summary], method = "spearman"),
    biobank_support, mean(sign(bio_all$pooled_effect[bio_valid]) == sign(meta_fit$pooled_effect[bio_valid])), cor(bio_all$pooled_effect[bio_valid], meta_fit$pooled_effect[bio_valid]), sum(bio_all$BH_FDR < 0.05, na.rm = TRUE),
    bio_sens_summary, g259_sens_summary, g259_labor_summary,
    max(tau_diff), max(effect_diff)
  ),
  stringsAsFactors = FALSE
)
write_csv(summary, file.path(interim, "phase1b_summary.csv"))

session <- sub("[ \\t]+$", "", capture.output(sessionInfo()))
writeLines(c(
  "Random seed: 20260809", "Primary meta method: vectorized intercept-only REML; modified Knapp-Hartung t inference", "REML validation: metafor 5.0.1 first 40 alphabetical eligible genes",
  paste0("Mapping version: HGNC_2026-08-09_faaeb6ae1e2a596b"), session
), file.path(qc_dir, "phase1b_session_info.txt"), useBytes = TRUE)

message("Phase 1B analysis complete. STABLE genes: ", stable_n)
