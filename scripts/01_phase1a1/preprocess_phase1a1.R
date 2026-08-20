#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(preprocessCore)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
raw_dir <- file.path(root, "data", "raw", "phase1a1")
interim <- file.path(root, "data", "interim", "phase1a1")
recon_dir <- file.path(interim, "reconstructed")
extract_dir <- file.path(interim, "raw_extract")
dir.create(recon_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

write_matrix <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(cbind(feature_id = rownames(x), x), con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

normalize_label <- function(x) gsub("[^a-z0-9]", "", tolower(x))
crosswalk <- read.csv(file.path(interim, "priority_sample_crosswalk.csv"), check.names = FALSE, stringsAsFactors = FALSE)

diagnostics <- list()

# GSE30186: the RAW tar contains platform manifests only. Reconstruct from the
# separate non-normalized BeadStudio AVG_Signal export using a control-free
# normexp background model (offset 16), quantile normalization, then log2.
g301 <- read.delim(gzfile(file.path(raw_dir, "GSE30186_non_normalized.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
avg_cols <- grep("\\.AVG_Signal$", colnames(g301), value = TRUE)
raw301 <- as.matrix(g301[, avg_cols, drop = FALSE])
storage.mode(raw301) <- "double"
rownames(raw301) <- g301$ID_REF
labels301 <- sub("\\.AVG_Signal$", "", avg_cols)
cw301 <- crosswalk[crosswalk$dataset == "GSE30186", , drop = FALSE]
match301 <- match(normalize_label(labels301), normalize_label(cw301$sample_title))
stopifnot(!anyNA(match301), ncol(raw301) == 12L)
colnames(raw301) <- cw301[["GSM/sample ID"]][match301]
bc301 <- backgroundCorrect(raw301, method = "normexp", normexp.method = "saddle", offset = 16)
qn301 <- normalizeBetweenArrays(bc301, method = "quantile")
recon301 <- log2(qn301)
write_matrix(recon301, file.path(recon_dir, "GSE30186_normexp16_quantile_log2_feature_level.tsv.gz"))
diagnostics[[length(diagnostics) + 1L]] <- data.frame(
  dataset = "GSE30186", submitted_feature_n = nrow(raw301), sample_n = ncol(raw301),
  raw_min = min(raw301, na.rm = TRUE), raw_max = max(raw301, na.rm = TRUE),
  raw_nonpositive_n = sum(raw301 <= 0, na.rm = TRUE), control_definition_n = 887,
  control_intensity_n = 0, flag_above_4096_n = NA_integer_,
  reconstruction = "limma normexp saddle offset=16; quantile normalization; log2",
  stringsAsFactors = FALSE
)

# GSE10588: retain the 32,878 regular GPL2986 probes, reconstruct the submitted
# ABarray-style path from raw Signal by quantile normalization then log2.
g105_extract <- file.path(extract_dir, "GSE10588")
dir.create(g105_extract, recursive = TRUE, showWarnings = FALSE)
untar(file.path(raw_dir, "GSE10588_RAW.tar"), exdir = g105_extract)
files105 <- sort(list.files(g105_extract, pattern = "^GSM[0-9]+\\.txt\\.gz$", full.names = TRUE))
ann105 <- read.delim(gzfile(file.path(interim, "GPL2986_annotation.tsv.gz")), check.names = FALSE, stringsAsFactors = FALSE)
probe105 <- as.character(ann105$ID)
raw105 <- matrix(NA_real_, nrow = length(probe105), ncol = length(files105), dimnames = list(probe105, sub("\\.txt\\.gz$", "", basename(files105))))
flags105 <- raw105
control_rows105 <- integer(length(files105))
for (i in seq_along(files105)) {
  tab <- read.delim(gzfile(files105[[i]]), check.names = FALSE, stringsAsFactors = FALSE)
  idx <- match(probe105, as.character(tab$Probe_ID))
  stopifnot(!anyNA(idx))
  raw105[, i] <- as.numeric(tab$Signal[idx])
  flags105[, i] <- as.numeric(tab$Flags[idx])
  control_rows105[[i]] <- sum(!as.character(tab$Probe_ID) %in% probe105)
}
qn105 <- normalize.quantiles(raw105, keep.names = TRUE)
recon105 <- log2(qn105)
write_matrix(recon105, file.path(recon_dir, "GSE10588_raw_signal_quantile_log2_feature_level.tsv.gz"))
diagnostics[[length(diagnostics) + 1L]] <- data.frame(
  dataset = "GSE10588", submitted_feature_n = nrow(raw105), sample_n = ncol(raw105),
  raw_min = min(raw105, na.rm = TRUE), raw_max = max(raw105, na.rm = TRUE),
  raw_nonpositive_n = sum(raw105 <= 0, na.rm = TRUE), control_definition_n = max(control_rows105),
  control_intensity_n = max(control_rows105), flag_above_4096_n = sum(flags105 > 4096, na.rm = TRUE),
  reconstruction = "regular-probe raw Signal; quantile normalization on natural scale; log2; flag variants audited separately",
  stringsAsFactors = FALSE
)

# GSE43942: use the 135,096 archived BLOCK1 PM probes. The design has three
# probes for 45,031 transcripts, two for one transcript and one for one
# transcript; retain the incomplete probe sets rather than silently dropping
# submitted transcript targets.
# RMA background correction, natural-scale quantile normalization, log2 and
# median-polish transcript summarization.
g439_extract <- file.path(extract_dir, "GSE43942")
dir.create(g439_extract, recursive = TRUE, showWarnings = FALSE)
untar(file.path(raw_dir, "GSE43942_RAW.tar"), exdir = g439_extract)
files439 <- sort(list.files(g439_extract, pattern = "^GSM[0-9]+.*_532\\.pair\\.gz$", full.names = TRUE))
first439 <- read.delim(gzfile(files439[[1]]), skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
keep439 <- first439$GENE_EXPR_OPTION == "BLOCK1"
probe439 <- as.character(first439$PROBE_ID[keep439])
seq439 <- as.character(first439$SEQ_ID[keep439])
raw439 <- matrix(NA_real_, nrow = length(probe439), ncol = length(files439), dimnames = list(probe439, sub("_.*$", "", basename(files439))))
raw439[, 1] <- as.numeric(first439$PM[keep439])
for (i in seq_along(files439)[-1]) {
  tab <- read.delim(gzfile(files439[[i]]), skip = 1, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- tab$GENE_EXPR_OPTION == "BLOCK1"
  stopifnot(identical(probe439, as.character(tab$PROBE_ID[keep])), identical(seq439, as.character(tab$SEQ_ID[keep])))
  raw439[, i] <- as.numeric(tab$PM[keep])
}
probe_count439 <- table(seq439)
stopifnot(length(unique(seq439)) == 45033L, ncol(raw439) == 12L,
          sum(probe_count439 == 3L) == 45031L,
          sum(probe_count439 == 2L) == 1L,
          sum(probe_count439 == 1L) == 1L)
bc439 <- rma.background.correct(raw439)
colnames(bc439) <- colnames(raw439)
qn439 <- normalize.quantiles(bc439, keep.names = TRUE)
colnames(qn439) <- colnames(raw439)
log439 <- log2(qn439)
group439 <- as.integer(factor(seq439, levels = unique(seq439)))
fits439 <- subrcModelMedianPolish(log439, group439)
recon439 <- matrix(NA_real_, nrow = length(fits439), ncol = ncol(log439), dimnames = list(unique(seq439), colnames(log439)))
for (i in seq_along(fits439)) recon439[i, ] <- fits439[[i]]$Estimates[seq_len(ncol(log439))]
write_matrix(recon439, file.path(recon_dir, "GSE43942_pair_rma_quantile_medianpolish_log2_transcript_level.tsv.gz"))
diagnostics[[length(diagnostics) + 1L]] <- data.frame(
  dataset = "GSE43942", submitted_feature_n = nrow(recon439), sample_n = ncol(recon439),
  raw_min = min(raw439, na.rm = TRUE), raw_max = max(raw439, na.rm = TRUE),
  raw_nonpositive_n = sum(raw439 <= 0, na.rm = TRUE), control_definition_n = nrow(first439) - sum(keep439),
  control_intensity_n = nrow(first439) - sum(keep439), flag_above_4096_n = NA_integer_,
  reconstruction = "PAIR PM; preprocessCore RMA background; natural-scale quantile; log2; median-polish by transcript (45031x3 probes, 1x2 probes, 1x1 probe)",
  stringsAsFactors = FALSE
)

write.csv(do.call(rbind, diagnostics), file.path(interim, "reconstruction_diagnostics.csv"), row.names = FALSE, na = "UNRESOLVED")
cat("Phase 1A.1 R reconstruction complete\n")
