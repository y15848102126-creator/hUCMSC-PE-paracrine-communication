#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args))normalizePath(args[1],winslash="/",mustWork=TRUE)else normalizePath(".",winslash="/");setwd(root)
out<-"results/03_phase3/figures";dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme<-function(){par(mar=c(5,5,2,1),las=1,bty="l",family="sans",cex.axis=.9,cex.lab=1)}
png_open<-function(name,w=1400,h=1000){png(file.path(out,name),width=w,height=h,res=160,bg="white");theme()}

r<-fread("results/03_phase3/baseline/baseline_sender_robustness.csv")
c<-fread("results/03_phase3/baseline/cross_dataset_sender_concordance.csv")
l<-fread("results/03_phase3/licensing/licensing_ligand_classification.csv")
h<-fread("results/03_phase3/sender/sender_evidence_hierarchy.csv")

png_open("phase3_baseline_classification.png")
par(mar=c(10,5,2,1))
ord<-c("ROBUST_BASELINE_SENDER","MULTIDATASET_LOW_EXPRESSION","DATASET_SPECIFIC","DONOR_VARIABLE","NOT_RELIABLY_EXPRESSED")
n<-table(r$baseline_classification);vals<-setNames(rep(0,length(ord)),ord);vals[names(n)]<-as.numeric(n)
barplot(vals,col=c("#176B87","#64CCC5","#DAA520","#B96D40","#B8BDC5"),ylab="Frozen ligand-universe genes",cex.names=.72,las=2)
title("Cross-dataset baseline sender classification",adj=0)
dev.off()

png_open("phase3_cross_dataset_rank_concordance.png")
cols<-ifelse(r$baseline_classification=="ROBUST_BASELINE_SENDER","#176B87AA","#9DA4AE55")
plot(r$gse182_median_rank,r$gse199_median_rank,pch=16,col=cols,xlim=c(0,1),ylim=c(0,1),xlab="GSE182158 UC donor median rank",ylab="GSE199071 donor median rank")
abline(0,1,lty=2,col="#555555");legend("topleft",legend=c("ROBUST_BASELINE_SENDER","Other"),pch=16,col=c("#176B87AA","#9DA4AE88"),bty="n")
title(sprintf("Independent baseline rank concordance (Spearman rho = %.3f)",unique(c$overall_spearman_rho)),adj=0)
dev.off()

png_open("phase3_licensing_classification.png")
par(mar=c(10,5,2,1))
ord2<-c("LICENSING_CONSISTENT_UP","LICENSING_CONSISTENT_DOWN","DONOR_DEPENDENT","PASSAGE_DEPENDENT","NO_CLEAR_LICENSING_EFFECT")
n2<-table(l$licensing_classification);v2<-setNames(rep(0,length(ord2)),ord2);v2[names(n2)]<-as.numeric(n2)
barplot(v2,col=c("#B33A3A","#2F6FB0","#D29F32","#8A5FA8","#B8BDC5"),ylab="Frozen ligand-universe genes",cex.names=.68,las=2)
title("Licensing classification (two donors; three valid strata)",adj=0)
dev.off()

png_open("phase3_donor_licensing_concordance.png")
pal<-c(LICENSING_CONSISTENT_UP="#B33A3AAA",LICENSING_CONSISTENT_DOWN="#2F6FB0AA",DONOR_DEPENDENT="#D29F3288",PASSAGE_DEPENDENT="#8A5FA866",NO_CLEAR_LICENSING_EFFECT="#AAB0B955")
plot(l$donor1_effect,l$donor2_effect,pch=16,col=pal[l$licensing_classification],xlim=range(l$donor1_effect,finite=TRUE),ylim=range(l$donor2_effect,finite=TRUE),xlab="Donor 1 P5 licensing log2 effect",ylab="Donor 2 median P2/P5 licensing log2 effect")
abline(h=0,v=0,lty=2,col="#777777");legend("bottomright",legend=names(pal),pch=16,col=pal,bty="n",cex=.72)
title("Donor-level inflammatory-licensing concordance",adj=0)
dev.off()

png_open("phase3_sender_evidence_levels.png")
ord3<-c("S1","S2","S3","S4");n3<-table(h$sender_evidence_level);v3<-setNames(rep(0,4),ord3);v3[names(n3)]<-as.numeric(n3)
barplot(v3,col=c("#0B6E4F","#368F8B","#F0A202","#B8BDC5"),ylab="Frozen ligand-universe genes",xlab="Independent sender evidence level")
title("Phase 3 sender evidence hierarchy",adj=0)
dev.off()

si<-sub("[[:space:]]+$","",capture.output(sessionInfo()));writeLines(si,file.path(out,"figure_session_info.txt"))
message("Phase 3 figures complete")
