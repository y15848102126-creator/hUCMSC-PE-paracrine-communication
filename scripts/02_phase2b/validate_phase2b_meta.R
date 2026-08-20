#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table);library(metafor)})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
setwd(root)

models <- fread("results/02_phase2b/scores/cohort_program_scores_summary.csv")
meta <- fread("results/02_phase2b/meta/program_gene_set_meta_analysis.csv")
diff_mu <- numeric(); diff_tau <- numeric()
for (id in meta$hypothesis_id) {
  d <- models[hypothesis_id == id]
  fit <- rma.uni(yi = d$beta_disease, sei = d$SE, method = "REML", test = "knha")
  saved <- meta[hypothesis_id == id]
  diff_mu <- c(diff_mu, abs(as.numeric(fit$b) - saved$summary_standardized_beta))
  diff_tau <- c(diff_tau, abs(fit$tau2 - saved$tau2))
}
cat(sprintf("Phase 2B REML cross-check differences: max beta %.6g; max tau2 %.6g\n", max(diff_mu), max(diff_tau)))
stopifnot(max(diff_mu) < 1e-5, max(diff_tau) < 1e-5)
cat("Phase 2B REML cross-check: PASS\n")
