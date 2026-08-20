#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
lib <- file.path(root, "data", "interim", "phase1a1", "Rlib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib)
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
bioc <- c("limma")
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, lib = lib)
if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor", repos = "https://cloud.r-project.org", lib = lib)
stopifnot(requireNamespace("limma", quietly = TRUE), requireNamespace("metafor", quietly = TRUE))
cat("limma=", as.character(packageVersion("limma")), "\n", sep = "")
cat("metafor=", as.character(packageVersion("metafor")), "\n", sep = "")
