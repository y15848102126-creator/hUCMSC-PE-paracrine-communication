# DEPRECATED / HISTORICAL ONLY: Phase 2A edgeR count-likelihood receiver inference is not authoritative. Use scripts/02_phase2a2/run_phase2a2_corrected_analysis.R.
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
script_arg <- commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))][1]
root <- normalizePath(file.path(dirname(sub("^--file=","",script_arg)),"../.."),winslash="/",mustWork=TRUE)
setwd(root)
figdir <- "results/02_phase2a/figures"; dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
patients <- fread("results/02_phase2a/metadata/patient_registry.csv")
pbreg <- fread("results/02_phase2a/pseudobulk/pseudobulk_registry.csv")
mds <- fread("results/02_phase2a/qc/pseudobulk_mds_coordinates.csv")
programs <- fread("results/02_phase2a/programs/pe_cellstate_programs.csv")
de <- fread("results/02_phase2a/DE/celltype_DE_summary.csv")

png(file.path(figdir,"A_cells_per_patient.png"),1200,760,res=130); par(mar=c(10,5,4,2))
cols<-c(CONTROL="#4C78A8",PE="#E45756")[patients$disease_status]
barplot(patients$cell_count,names.arg=patients$patient_id,las=2,col=cols,main="Published cells per patient",ylab="Cells",cex.names=.65);legend("topright",legend=c("Control","PE"),fill=c("#4C78A8","#E45756"),bty="n");dev.off()

heat<-dcast(pbreg,patient_id~harmonized_annotation,value.var="cell_count");hm<-as.matrix(heat[,-1]);rownames(hm)<-heat$patient_id
png(file.path(figdir,"B_cells_patient_celltype_heatmap.png"),1500,1150,res=130);par(mar=c(13,9,4,2))
image(t(log10(hm+1)[nrow(hm):1,]),axes=FALSE,col=hcl.colors(40,"YlOrRd"),main="log10(cells + 1) per patient × cell type");axis(1,at=seq(0,1,length.out=ncol(hm)),labels=gsub("IMMUNE_PROGENITOR_OR_PROLIFERATING","IMMUNE_PROGENITOR",colnames(hm)),las=2,cex.axis=.65);axis(2,at=seq(0,1,length.out=nrow(hm)),labels=rev(rownames(hm)),las=2,cex.axis=.65);dev.off()

png(file.path(figdir,"C_pseudobulk_library_sizes.png"),1250,820,res=130);par(mar=c(12,6,4,2))
pbreg[,plot_annotation:=fcase(harmonized_annotation=="IMMUNE_PROGENITOR_OR_PROLIFERATING","IMMUNE_PROGENITOR",harmonized_annotation=="VASCULAR_SMOOTH_MUSCLE","VASCULAR_SM",default=harmonized_annotation)]
boxplot(log10(pseudobulk_library_umi+1)~plot_annotation,data=pbreg,las=2,col="#72B7B2",ylab="log10(count-matrix library UMI + 1)",main="Patient pseudobulk library sizes",cex.axis=.7,xlab="");abline(h=log10(10000+1),lty=2,col="red");dev.off()

top_ct<-de[contrast%in%c("EOPE","LOPE"),.(n=sum(n_PE+n_control)),by=celltype][order(-n)][1:min(6,.N),celltype]
zall<-mds[contrast%in%c("EOPE","LOPE")&celltype%in%top_ct]
png(file.path(figdir,"D_patient_pseudobulk_MDS.png"),1500,950,res=130);par(mfrow=c(2,3),mar=c(4,4,3,1))
for(ct in top_ct){z<-zall[celltype==ct];plot(z$MDS1,z$MDS2,col=ifelse(z$disease=="PE","#E45756","#4C78A8"),pch=19,xlab="MDS1",ylab="MDS2",main=ct);text(z$MDS1,z$MDS2,labels=z$patient_id,pos=3,cex=.42)};dev.off()

pc<-unique(programs[contrast%in%c("EOPE","LOPE")&BH_FDR<.05,.(celltype,pathway,classification)])[, .N,by=.(celltype,classification)]
png(file.path(figdir,"E_significant_program_counts.png"),1400,900,res=130);par(mar=c(12,6,4,10),xpd=NA)
widepc<-dcast(pc,celltype~classification,value.var="N",fill=0);m<-as.matrix(widepc[,-1]);rownames(m)<-widepc$celltype;pal<-hcl.colors(ncol(m),"Set2");barplot(t(m),las=2,col=pal,main="Significant subtype-stratified program classifications",ylab="Unique gene sets",cex.names=.7);legend("topright",inset=c(-.22,0),legend=colnames(m),fill=pal,cex=.72,bty="n");dev.off()

top<-head(unique(programs[classification=="SHARED_PE"&contrast%in%c("EOPE","LOPE"),.(celltype,pathway,gene_set)]),30)
png(file.path(figdir,"F_shared_program_heatmap.png"),1600,1100,res=130);par(mar=c(7,18,4,2))
if(nrow(top)){z<-merge(programs[contrast%in%c("EOPE","LOPE"),.(celltype,pathway,contrast,NES)],top,by=c("celltype","pathway"));z[,short:=substr(gsub("^(HALLMARK_|REACTOME_|GOBP_)","",gene_set),1,37)];z[,cell_short:=fifelse(celltype=="PLACENTAL_STROMAL","P_STROMAL",celltype)];z[,label:=paste(cell_short,short,sep=" | ")];w<-dcast(z,label~contrast,value.var="NES");mm<-as.matrix(w[,-1]);rownames(mm)<-w$label;image(t(mm[nrow(mm):1,,drop=FALSE]),axes=FALSE,col=hcl.colors(51,"Blue-Red 3",rev=TRUE),main="Shared PE programs: subtype signed score");axis(1,at=seq(0,1,length.out=ncol(mm)),labels=colnames(mm),las=2);axis(2,at=seq(0,1,length.out=nrow(mm)),labels=rev(rownames(mm)),las=2,cex.axis=.58)}else{plot.new();text(.5,.5,"No SHARED_PE programs under frozen criteria")};dev.off()
