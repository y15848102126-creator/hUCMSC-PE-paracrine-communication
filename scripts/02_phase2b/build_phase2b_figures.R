#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
setwd(root)

out <- "results/02_phase2b/figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("GSE75010_BIOBANK", "GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942")
mods <- sprintf("PROGRAM_MODULE_%02d", c(1:7, 9:11))
short_mod <- sub("PROGRAM_MODULE_", "M", mods)

models <- fread("results/02_phase2b/scores/cohort_program_scores_summary.csv")
meta <- fread("results/02_phase2b/meta/program_gene_set_meta_analysis.csv")
mv <- fread("results/02_phase2b/meta/program_module_validation.csv")
loco <- fread("results/02_phase2b/robustness/program_leave_one_cohort_out.csv")
ev <- fread("results/02_phase2b/evidence/updated_receiver_evidence_hierarchy.csv")

safe_png <- function(path, width = 1800, height = 1200) {
  png(path, width = width, height = height, res = 170, type = "cairo-png")
}
heat_cols <- colorRampPalette(c("#2B6CB0", "#F7FAFC", "#C53030"))(101)

# A. Each tile is the fraction of a module's frozen constituents agreeing with
# the frozen scRNA direction in that cohort. No pooled module gene union is made.
aligned <- copy(models)
aligned[, expected_agreement := as.integer(agrees_scRNA_direction == "YES")]
mh <- aligned[, .(agreement = mean(expected_agreement)), by = .(program_module, cohort)]
mat <- matrix(NA_real_, nrow = length(mods), ncol = length(cohorts), dimnames = list(short_mod, cohorts))
for (i in seq_len(nrow(mh))) mat[match(mh$program_module[i], mods), match(mh$cohort[i], cohorts)] <- mh$agreement[i]
safe_png(file.path(out, "phase2b_module_by_cohort_direction_heatmap.png"), 1900, 1200)
par(mar = c(9, 7, 4, 2))
image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), zlim = c(0, 1),
      col = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(101), axes = FALSE,
      xlab = "", ylab = "", main = "Frozen module direction agreement by independent bulk cohort")
axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = .78)
axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 1)
for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
  z <- mat[nrow(mat) - i + 1, j]
  text(j, i, sprintf("%.0f%%", 100 * z), cex = .75, col = ifelse(z < .18 || z > .82, "white", "#1A202C"))
}
mtext("Red = 0%, white = 50%, blue = 100% expected direction. Bulk placenta does not validate cell-type localization.", side = 1, line = 7.4, cex = .72)
dev.off()

# B. Representative is frozen outcome-independently as the lexicographically
# first constituent name within each module.
rep_ids <- mv$representative_hypothesis_id
fm <- meta[match(rep_ids, hypothesis_id)]
ord <- rev(seq_len(nrow(fm)))
safe_png(file.path(out, "phase2b_representative_constituent_forest.png"), 1900, 1250)
par(mar = c(7, 14, 4, 3))
xlim <- range(c(fm$CI_lower, fm$CI_upper, 0), finite = TRUE)
plot(NA, xlim = xlim, ylim = c(.5, nrow(fm) + .5), yaxt = "n", xlab = "Random-effects summary standardized disease beta", ylab = "",
     main = "Outcome-independent representative constituent per frozen module")
abline(v = 0, col = "#718096", lty = 2)
segments(fm$CI_lower[ord], seq_along(ord), fm$CI_upper[ord], seq_along(ord), col = "#4A5568", lwd = 2)
cols <- ifelse(fm$classification[ord] == "BULK_DIRECTIONAL_SUPPORT", "#2B6CB0", "#A0AEC0")
points(fm$summary_standardized_beta[ord], seq_along(ord), pch = 19, col = cols, cex = 1.05)
axis(2, at = seq_along(ord), labels = sub("PROGRAM_MODULE_", "M", fm$program_module[ord]), las = 1)
legend("bottomright", legend = c("Directional support", "Not supported"), col = c("#2B6CB0", "#A0AEC0"), pch = 19, bty = "n")
mtext("Representative rule: lexicographically first frozen constituent; selection did not use bulk outcomes.", side = 1, line = 5.2, cex = .72)
dev.off()

# C. All sets that reached the frozen directional/robust/heterogeneous support
# family are shown; none reached robust in this analysis.
potential <- meta[classification %in% c("BULK_ROBUST_SUPPORT", "BULK_DIRECTIONAL_SUPPORT", "BULK_HETEROGENEOUS_SUPPORT"), hypothesis_id]
lp <- loco[hypothesis_id %in% potential]
if (!"scRNA_direction" %in% names(lp)) lp <- merge(lp, meta[, .(hypothesis_id, scRNA_direction)], by = "hypothesis_id", all.x = TRUE)
lp[, aligned_beta := summary_standardized_beta * ifelse(scRNA_direction == "UP_IN_PE_PROGRAM", 1, -1)]
lmat <- matrix(NA_real_, nrow = length(potential), ncol = length(cohorts), dimnames = list(potential, cohorts))
for (i in seq_len(nrow(lp))) lmat[match(lp$hypothesis_id[i], potential), match(lp$omitted_cohort[i], cohorts)] <- lp$aligned_beta[i]
lim <- max(abs(lmat), na.rm = TRUE)
safe_png(file.path(out, "phase2b_leave_one_cohort_out_robustness.png"), 1900, 1200)
par(mar = c(9, 6, 4, 2))
image(seq_len(ncol(lmat)), seq_len(nrow(lmat)), t(lmat[nrow(lmat):1, , drop = FALSE]), zlim = c(-lim, lim), col = heat_cols,
      axes = FALSE, xlab = "Omitted cohort", ylab = "", main = "LOCO expected-direction-aligned effects for directionally supported sets")
axis(1, at = seq_len(ncol(lmat)), labels = colnames(lmat), las = 2, cex.axis = .78)
axis(2, at = seq_len(nrow(lmat)), labels = rev(rownames(lmat)), las = 1, cex.axis = .85)
abline(h = seq(.5, nrow(lmat) + .5, 1), v = seq(.5, ncol(lmat) + .5, 1), col = "white", lwd = .5)
mtext("Red = expected direction; blue = opposite. Positive aligned beta retains the frozen scRNA direction after omission.", side = 1, line = 7.4, cex = .72)
dev.off()

# D. Evidence domains remain separate. Module 08 is included only to display its
# frozen hold status and was not tested in Phase 2B.
ev <- ev[match(sprintf("PROGRAM_MODULE_%02d", 1:11), program_module)]
bulk_code <- fcase(ev$INDEPENDENT_BULK_PROGRAM_SUPPORT == "BULK_MODULE_SUPPORTED", 2,
                   ev$INDEPENDENT_BULK_PROGRAM_SUPPORT == "BULK_MODULE_DIRECTIONAL", 1,
                   ev$INDEPENDENT_BULK_PROGRAM_SUPPORT == "NOT_TESTED_HOLD_EXTERNAL_DISCORDANCE", -1,
                   default = 0)
ext_code <- fcase(ev$EXTERNAL_SCRNA_SUPPORT == "YANG_LOPE_DIRECTIONAL_SUPPORT", 2,
                  ev$EXTERNAL_SCRNA_SUPPORT == "DIRECTIONALLY_DISCORDANT", -1, default = 0)
adm_code <- rep(2, nrow(ev))
emat <- cbind(adm_code, ext_code, bulk_code)
rownames(emat) <- sub("PROGRAM_MODULE_", "M", ev$program_module)
colnames(emat) <- c("Corrected Admati", "External scRNA", "Independent bulk")
ecol <- c("#B2182B", "#E2E8F0", "#90CDF4", "#2166AC"); names(ecol) <- c("-1", "0", "1", "2")
safe_png(file.path(out, "phase2b_scrna_bulk_evidence_matrix.png"), 1900, 1250)
par(mar = c(10, 7, 4, 13))
plot(c(.5, ncol(emat) + .5), c(.5, nrow(emat) + .5), type = "n", axes = FALSE, xlab = "", ylab = "",
     main = "Receiver evidence domains kept separate")
for (i in seq_len(nrow(emat))) for (j in seq_len(ncol(emat))) rect(j - .48, nrow(emat) - i + .52, j + .48, nrow(emat) - i + 1.48, col = ecol[as.character(emat[i, j])], border = "white")
axis(1, at = seq_len(ncol(emat)), labels = colnames(emat), las = 2)
axis(2, at = seq_len(nrow(emat)), labels = rev(rownames(emat)), las = 1)
legend("right", inset = c(-.27, 0), xpd = TRUE,
       legend = c("Strong/corrected support", "Directional support", "Not evaluable/not supported", "Discordant/held"),
       fill = c("#2166AC", "#90CDF4", "#E2E8F0", "#B2182B"), bty = "n", cex = .72)
mtext("Bulk evidence concerns placental tissue-level program direction, never cellular origin.", side = 1, line = 8.3, cex = .72)
dev.off()

message("Phase 2B analytical figure previews written")
