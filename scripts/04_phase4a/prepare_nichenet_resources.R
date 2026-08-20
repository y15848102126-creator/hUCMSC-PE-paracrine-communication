#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
out <- file.path(root,"data/interim/phase4a")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

candidates <- read.csv(file.path(root,"results/03_phase3/sender/frozen_phase4_sender_candidates.csv"),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(candidates)==214L)

lt <- readRDS(file.path(root,"data/raw/phase4a/ligand_target_matrix_nsga2r_final.rds"))
keep <- intersect(candidates$gene,colnames(lt))
stopifnot(length(keep)>=200L)
sub <- as.matrix(lt[,keep,drop=FALSE])
writeLines(rownames(sub),file.path(out,"nichenet_target_rows.txt"),useBytes=TRUE)
writeLines(colnames(sub),file.path(out,"nichenet_target_columns.txt"),useBytes=TRUE)
con <- file(file.path(out,"nichenet_target_subset_float64.bin"),"wb")
writeBin(as.double(sub),con,size=8,endian="little")
close(con)
writeLines(c(sprintf("rows=%d",nrow(sub)),sprintf("columns=%d",ncol(sub)),"storage=float64_column_major_little_endian"),file.path(out,"nichenet_target_subset_manifest.txt"))

lr <- as.data.frame(readRDS(file.path(root,"data/raw/phase3/lr_network_human_21122021.rds")),stringsAsFactors=FALSE)
write.table(lr,file.path(out,"nichenet_lr_network.tsv"),sep="\t",row.names=FALSE,col.names=TRUE,quote=TRUE,na="")
message(sprintf("Prepared NicheNet resources: %d targets x %d candidate ligands; %d LR edges",nrow(sub),ncol(sub),nrow(lr)))
