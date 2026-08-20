#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
out <- file.path(root,"results/04_phase4a/figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

ev <- read.csv(file.path(root,"results/04_phase4a/integration/sender_receiver_evidence_matrix.csv"),check.names=FALSE,stringsAsFactors=FALSE)
hier <- read.csv(file.path(root,"results/04_phase4a/integration/phase4a_candidate_hierarchy.csv"),check.names=FALSE,stringsAsFactors=FALSE)

module_order <- sprintf("PROGRAM_MODULE_%02d",c(1,4,5,6,7,10))
module_label <- paste0("M",c(1,4,5,6,7,10))
compatible <- sapply(module_order,function(x) sum(ev$program_module==x & ev$competent_receptor_n>0))
target <- sapply(module_order,function(x) sum(ev$program_module==x & ev$target_compatibility_class=="TARGET_COMPATIBLE"))
reversal <- sapply(module_order,function(x) sum(ev$program_module==x & ev$signed_reversal_class=="REVERSAL_SUPPORTED"))
mat <- rbind(Receptor_competent=compatible,Target_compatible=target,Signed_reversal=reversal)
png(file.path(out,"phase4a_module_evidence_counts.png"),width=1400,height=850,res=150)
par(mar=c(5,5,2,1),las=1)
barplot(mat,beside=TRUE,names.arg=module_label,col=c("#4477AA","#66CCEE","#228833"),ylab="Frozen sender candidates (n)",main="Blinded Phase 4A evidence counts by receiver module")
legend("topright",legend=c("Competent receptor","Target compatible (unsigned)","Signed reversal"),fill=c("#4477AA","#66CCEE","#228833"),bty="n")
mtext("Counts are computational evidence categories, not therapeutic proof",side=1,line=3.5,cex=.8)
dev.off()

classes <- c("REVERSAL_SUPPORTED","DISEASE_CONCORDANT_POTENTIAL","SIGNED_EVIDENCE_INSUFFICIENT")
signed <- sapply(module_order,function(x) table(factor(ev$signed_reversal_class[ev$program_module==x],levels=classes)))
png(file.path(out,"phase4a_signed_evidence_landscape.png"),width=1400,height=850,res=150)
par(mar=c(5,5,2,1),las=1)
barplot(signed,names.arg=module_label,space=0.55,cex.names=0.78,col=c("#228833","#CC6677","#BBBBBB"),ylab="Ligand-receiver axes (n)",main="Signed evidence remains separate from unsigned compatibility")
legend("topright",legend=c("Reversal supported","Disease-concordant potential","Signed evidence insufficient"),fill=c("#228833","#CC6677","#BBBBBB"),bty="n")
dev.off()

tier_levels <- c("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE","TIER_B_COMPATIBILITY_CANDIDATE","TIER_C_EXTENDED_EXTRACELLULAR","NOT_PRIORITIZED")
tier_short <- c("Tier A","Tier B","Tier C","Not prioritized")
by_sender <- table(factor(hier$best_phase4a_tier,levels=tier_levels),factor(hier$sender_evidence_level,levels=c("S1","S2")))
png(file.path(out,"phase4a_sender_class_descriptive.png"),width=1300,height=800,res=150)
par(mar=c(6,5,2,1),las=1)
barplot(by_sender,beside=TRUE,names.arg=c("S1","S2"),col=c("#332288","#88CCEE","#DDCC77","#BBBBBB"),ylab="Frozen sender candidates (n)",main="Descriptive Phase 4A hierarchy by frozen sender class")
legend("topright",legend=tier_short,fill=c("#332288","#88CCEE","#DDCC77","#BBBBBB"),bty="n")
mtext("S1 is not assigned priority over S2",side=1,line=4.5,cex=.8)
dev.off()

cat("PHASE4A_FIGURES_OK\n")
