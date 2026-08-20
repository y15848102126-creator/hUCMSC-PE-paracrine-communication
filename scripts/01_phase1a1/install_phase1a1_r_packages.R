#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/")
lib <- file.path(root, "data", "interim", "phase1a1", "Rlib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib)
}
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
needed <- c("limma", "preprocessCore")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE, lib = lib)
stopifnot(all(vapply(needed, requireNamespace, logical(1), quietly = TRUE)))
cat(paste(needed, vapply(needed, function(x) as.character(packageVersion(x)), character(1)), sep = "="), sep = "\n")
