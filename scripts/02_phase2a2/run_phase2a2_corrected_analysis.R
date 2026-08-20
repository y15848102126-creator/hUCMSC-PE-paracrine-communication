#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[1], winslash="/", mustWork=TRUE) else normalizePath(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "../.."), winslash="/", mustWork=TRUE)
setwd(root)
set.seed(20260809)

out <- "results/02_phase2a2"
dir.create(file.path(out,"corrected_analysis"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"evidence"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"figures"),recursive=TRUE,showWarnings=FALSE)
write_csv <- function(x,p) fwrite(x,p,na="",quote=TRUE)

message("Reading sums of public normalized-and-ceiled values")
d <- fread("data/interim/phase2a/admati_harmonized_pseudobulk_counts.csv",check.names=FALSE)
gene_input <- d[[1]]
mat <- as.matrix(d[,-1]); storage.mode(mat) <- "double"
sums <- rowsum(mat,group=gene_input,reorder=FALSE)
rm(d,mat,gene_input); invisible(gc())

elig <- fread("results/02_phase2a/metadata/pseudobulk_eligibility.csv")
elig <- elig[contrast %in% c("EOPE","LOPE") & include_in_contrast=="YES"]
shared <- fread("results/02_phase2a/programs/shared_pe_programs.csv")
shared[,legacy_direction:=ifelse(NES_EOPE>0,"UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM")]
modules_legacy <- fread("results/02_phase2a1/shared_program_modules.csv")

read_gmt <- function(path,collection){
  lines <- readLines(path,warn=FALSE)
  ans <- lapply(lines,function(line){z<-strsplit(line,"\t",fixed=TRUE)[[1]];unique(z[-c(1,2)])})
  names(ans) <- paste0(collection,"::",vapply(lines,function(line)strsplit(line,"\t",fixed=TRUE)[[1]][1],character(1)))
  ans
}
resources <- "data/raw/phase2a_resources"
pathways <- c(
  read_gmt(file.path(resources,"h.all.v2026.1.Hs.symbols.gmt"),"HALLMARK"),
  read_gmt(file.path(resources,"c2.cp.reactome.v2026.1.Hs.symbols.gmt"),"REACTOME"),
  read_gmt(file.path(resources,"c5.go.bp.v2026.1.Hs.symbols.gmt"),"GOBP")
)
stopifnot(length(pathways)==9427L)

gene_rows <- list(); tier1_rows <- list(); tier2_rows <- list(); diagnostic_rows <- list()
models <- new.env(parent=emptyenv())
for(con in c("EOPE","LOPE")){
  for(ct in sort(unique(elig[contrast==con]$harmonized_annotation))){
    ee <- elig[contrast==con & harmonized_annotation==ct]
    ee[,matrix_column:=paste(patient_id,harmonized_annotation,sep="__")]
    stopifnot(all(ee$matrix_column %in% colnames(sums)))
    x <- sums[,ee$matrix_column,drop=FALSE]
    x <- sweep(x,2,ee$cell_count,"/")
    logexpr <- log2(x+1)
    keep <- rowSums(x>0)>=3L
    logexpr <- logexpr[keep,,drop=FALSE]
    group <- factor(ifelse(ee$group %in% c("EOPE","LOPE"),"PE","CONTROL"),levels=c("CONTROL","PE"))
    design <- model.matrix(~group)
    fit <- eBayes(lmFit(logexpr,design),trend=TRUE,robust=TRUE)
    coef_name <- "groupPE"
    tt <- topTable(fit,coef=coef_name,number=Inf,sort.by="none",adjust.method="BH")
    se <- fit$stdev.unscaled[,coef_name]*sqrt(fit$s2.post)
    gr <- data.table(
      gene=rownames(tt),contrast=con,celltype=ct,corrected_log2_mean_difference=tt$logFC,
      SE=se,moderated_t=tt$t,P=tt$P.Value,BH_FDR=tt$adj.P.Val,
      n_PE=sum(group=="PE"),n_control=sum(group=="CONTROL"),
      source_url="Figshare:23264102.v1|config:phase2a2"
    )
    gene_rows[[paste(con,ct,sep="|")]] <- gr
    stat <- setNames(tt$t,rownames(tt))
    models[[paste(con,ct,sep="|")]] <- list(stat=stat,n_PE=sum(group=="PE"),n_control=sum(group=="CONTROL"))
    diagnostic_rows[[paste(con,ct,sep="|")]] <- data.table(
      contrast=con,celltype=ct,patient_n=length(group),n_PE=sum(group=="PE"),n_control=sum(group=="CONTROL"),
      gene_input_n=nrow(sums),estimable_gene_n=sum(keep),design_rank=qr(design)$rank,design_columns=ncol(design),
      residual_df=min(fit$df.residual),median_residual_sigma=median(fit$sigma),
      eligibility_reused_exactly="YES",negative_binomial_likelihood="NO",model_status="PASS",
      source_url="config/phase2a2_analysis.json|results/02_phase2a/metadata/pseudobulk_eligibility.csv")

    s <- shared[celltype==ct]
    if(nrow(s)){
      idx <- lapply(s$pathway,function(p)intersect(pathways[[p]],names(stat))); names(idx)<-s$pathway
      valid <- lengths(idx)>=15 & lengths(idx)<=500
      if(any(valid)){
        cp <- as.data.table(cameraPR(stat,idx[valid],inter.gene.cor=0.01,sort=FALSE),keep.rownames="pathway")
        cp <- merge(cp,s[,.(celltype,collection,gene_set,pathway,legacy_direction)],by="pathway",all.x=TRUE)
        cp[,`:=`(contrast=con,available_gene_n=lengths(idx)[pathway],corrected_effect_statistic=vapply(pathway,function(p)mean(stat[idx[[p]]]),numeric(1)),
          corrected_direction=ifelse(Direction=="Up","UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM"),P_value=PValue,
          n_PE=sum(group=="PE"),n_control=sum(group=="CONTROL"),method="cameraPR on corrected limma moderated t",
          statistic_name="MEAN_MEMBER_MODERATED_T",source_url="MSigDB:2026.1.Hs|config/phase2a2_analysis.json")]
        tier1_rows[[paste(con,ct,sep="|")]] <- cp
      }
    }

    idx_all <- lapply(pathways,function(g)intersect(g,names(stat)))
    valid_all <- lengths(idx_all)>=15 & lengths(idx_all)<=500
    cp2 <- as.data.table(cameraPR(stat,idx_all[valid_all],inter.gene.cor=0.01,sort=FALSE),keep.rownames="pathway")
    cp2[,`:=`(contrast=con,celltype=ct,collection=sub("::.*$","",pathway),gene_set=sub("^[^:]+::","",pathway),
      available_gene_n=lengths(idx_all)[pathway],corrected_effect_statistic=vapply(pathway,function(p)mean(stat[idx_all[[p]]]),numeric(1)),
      corrected_direction=ifelse(Direction=="Up","UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM"),P_value=PValue,
      n_PE=sum(group=="PE"),n_control=sum(group=="CONTROL"),interpretation="EXPLORATORY_CORRECTED_REDISCOVERY",
      source_url="MSigDB:2026.1.Hs|config/phase2a2_analysis.json")]
    cp2[,BH_FDR:=p.adjust(P_value,"BH")]
    tier2_rows[[paste(con,ct,sep="|")]] <- cp2
    message(sprintf("%s %s: %d donors, %d genes, %d sets",con,ct,length(group),nrow(tt),nrow(cp2)))
  }
}

genes <- rbindlist(gene_rows,use.names=TRUE,fill=TRUE)
setorder(genes,celltype,contrast,P)
write_csv(genes,file.path(out,"corrected_analysis/corrected_gene_statistics.csv"))
write_csv(rbindlist(diagnostic_rows),file.path(out,"corrected_analysis/corrected_model_diagnostics.csv"))

t1long <- rbindlist(tier1_rows,use.names=TRUE,fill=TRUE)
t1long[,BH_FDR_frozen20:=p.adjust(P_value,"BH",n=20),by=contrast]
keep_t1 <- c("celltype","collection","gene_set","pathway","legacy_direction","contrast","available_gene_n","corrected_effect_statistic","statistic_name","corrected_direction","Direction","P_value","BH_FDR_frozen20","n_PE","n_control","method","source_url")
t1w <- dcast(t1long[,..keep_t1],celltype+collection+gene_set+pathway+legacy_direction~contrast,value.var=c("available_gene_n","corrected_effect_statistic","corrected_direction","P_value","BH_FDR_frozen20","n_PE","n_control"))
t1w[,direction_agrees_EOPE:=ifelse(corrected_direction_EOPE==legacy_direction,"YES","NO")]
t1w[,direction_agrees_LOPE:=ifelse(corrected_direction_LOPE==legacy_direction,"YES","NO")]
t1w[,classification:=fcase(
  direction_agrees_EOPE=="YES" & direction_agrees_LOPE=="YES" & BH_FDR_frozen20_EOPE<0.05 & BH_FDR_frozen20_LOPE<0.05,"CORRECTED_SHARED_SUPPORT",
  direction_agrees_EOPE=="YES" & BH_FDR_frozen20_EOPE<0.05 & !(direction_agrees_LOPE=="YES" & BH_FDR_frozen20_LOPE<0.05),"EOPE_ONLY_SUPPORT",
  direction_agrees_LOPE=="YES" & BH_FDR_frozen20_LOPE<0.05 & !(direction_agrees_EOPE=="YES" & BH_FDR_frozen20_EOPE<0.05),"LOPE_ONLY_SUPPORT",
  direction_agrees_EOPE=="YES" & direction_agrees_LOPE=="YES","DIRECTION_ONLY",
  default="NOT_SUPPORTED")]
t1w[,`:=`(historical_hypothesis_label="LEGACY_COUNT_MODEL_DISCOVERY",tier="TIER_1_FROZEN_HYPOTHESIS_RETEST",
  correction_family="BH_WITHIN_FROZEN_20_SEPARATELY_BY_SUBTYPE",source_url="results/02_phase2a/programs/shared_pe_programs.csv|config/phase2a2_analysis.json")]
setorder(t1w,celltype,pathway)
write_csv(t1w,file.path(out,"corrected_analysis/frozen20_corrected_retest.csv"))

t2 <- rbindlist(tier2_rows,use.names=TRUE,fill=TRUE)
t2w <- dcast(t2,celltype+collection+gene_set+pathway~contrast,value.var=c("available_gene_n","corrected_effect_statistic","corrected_direction","P_value","BH_FDR"))
t2w[,classification:=fcase(
  BH_FDR_EOPE<0.05 & BH_FDR_LOPE<0.05 & corrected_direction_EOPE==corrected_direction_LOPE,"CORRECTED_SHARED_PE",
  BH_FDR_EOPE<0.05 & (is.na(BH_FDR_LOPE)|BH_FDR_LOPE>=0.10),"CORRECTED_EOPE_ENRICHED",
  BH_FDR_LOPE<0.05 & (is.na(BH_FDR_EOPE)|BH_FDR_EOPE>=0.10),"CORRECTED_LOPE_ENRICHED",
  (pmin(BH_FDR_EOPE,BH_FDR_LOPE,na.rm=TRUE)<0.05 & corrected_direction_EOPE!=corrected_direction_LOPE) |
    (BH_FDR_EOPE<0.05 & BH_FDR_LOPE>=0.05 & BH_FDR_LOPE<0.10) |
    (BH_FDR_LOPE<0.05 & BH_FDR_EOPE>=0.05 & BH_FDR_EOPE<0.10),"CORRECTED_UNSTABLE",
  default="CORRECTED_NOT_SIGNIFICANT")]
t2w[,`:=`(interpretation="EXPLORATORY_CORRECTED_REDISCOVERY",source_url="MSigDB:2026.1.Hs|config/phase2a2_analysis.json")]
t2w[,sort_FDR:=pmin(BH_FDR_EOPE,BH_FDR_LOPE,na.rm=TRUE)]
setorder(t2w,classification,celltype,sort_FDR)
t2w[,sort_FDR:=NULL]
write_csv(t2w,file.path(out,"corrected_analysis/corrected_program_rediscovery.csv"))

orig <- modules_legacy[record_type=="ORIGINAL_GENE_SET",.(program_module,module_label,celltype,frozen_direction,collection,gene_set,pathway)]
orig <- merge(orig,t1w[,.(celltype,pathway,corrected_direction_EOPE,corrected_direction_LOPE,BH_FDR_frozen20_EOPE,BH_FDR_frozen20_LOPE,direction_agrees_EOPE,direction_agrees_LOPE,classification)],by=c("celltype","pathway"),all.x=TRUE)
module_summary <- orig[,{
  agree_both <- direction_agrees_EOPE=="YES" & direction_agrees_LOPE=="YES"
  any_sig <- any(BH_FDR_frozen20_EOPE<0.05 | BH_FDR_frozen20_LOPE<0.05,na.rm=TRUE)
  prop <- mean(agree_both,na.rm=TRUE)
  status <- if(prop>=0.5 && any_sig) "CORRECTED_ADMATI_SUPPORT" else if(prop>=0.5) "DIRECTIONAL_CORRECTED_SUPPORT_ONLY" else "NOT_SUPPORTED_AFTER_CORRECTION"
  list(record_type="PROGRAM_MODULE",constituent_gene_set_n=.N,direction_agree_both_n=sum(agree_both,na.rm=TRUE),direction_agree_both_proportion=prop,
       any_constituent_fdr05=ifelse(any_sig,"YES","NO"),corrected_module_status=status,
       constituent_pathways=paste(pathway,collapse=";"),source_url="results/02_phase2a1/shared_program_modules.csv|config/phase2a2_analysis.json")
},by=.(program_module,module_label,celltype,frozen_direction)]
orig[,`:=`(record_type="ORIGINAL_GENE_SET",constituent_gene_set_n=NA_integer_,direction_agree_both_n=NA_integer_,direction_agree_both_proportion=NA_real_,any_constituent_fdr05="",corrected_module_status="CONSTITUENT",constituent_pathways="",source_url="results/02_phase2a1/shared_program_modules.csv|results/02_phase2a/programs/shared_pe_programs.csv")]
mods <- rbindlist(list(orig,module_summary),use.names=TRUE,fill=TRUE)
setorder(mods,program_module,record_type,pathway)
write_csv(mods,file.path(out,"corrected_analysis/corrected_program_modules.csv"))

png(file.path(out,"figures/A_frozen20_corrected_direction.png"),1500,900,res=130)
plotdata <- t1w; vals <- rbind(plotdata$corrected_effect_statistic_EOPE,plotdata$corrected_effect_statistic_LOPE)
colnames(vals) <- paste0(plotdata$celltype,"\n",sub("^(HALLMARK|REACTOME|GOBP)_","",plotdata$gene_set)); rownames(vals)<-c("EOPE","LOPE")
z <- pmax(pmin(vals,quantile(abs(vals),.95,na.rm=TRUE)), -quantile(abs(vals),.95,na.rm=TRUE))
par(mar=c(14,5,4,2)); image(seq_len(ncol(z)),seq_len(nrow(z)),t(z),col=hcl.colors(101,"Blue-Red 3",rev=TRUE),axes=FALSE,xlab="",ylab="",main="Corrected mean member moderated-t: frozen 20"); axis(1,seq_len(ncol(z)),colnames(z),las=2,cex.axis=.55);axis(2,seq_len(nrow(z)),rownames(z),las=1);abline(v=seq(.5,ncol(z)+.5,1),col="white");dev.off()

png(file.path(out,"figures/B_tier1_classification.png"),1100,750,res=130)
tb <- table(factor(t1w$classification,levels=c("CORRECTED_SHARED_SUPPORT","EOPE_ONLY_SUPPORT","LOPE_ONLY_SUPPORT","DIRECTION_ONLY","NOT_SUPPORTED")))
par(mar=c(9,5,4,2));barplot(tb,las=2,col=hcl.colors(length(tb),"Set 2"),ylab="Frozen hypotheses",main="Corrected classification of 20 legacy hypotheses");dev.off()

png(file.path(out,"figures/C_corrected_module_status.png"),1200,750,res=130)
ms <- module_summary[order(program_module)]; cols <- c(CORRECTED_ADMATI_SUPPORT="#59A14F",DIRECTIONAL_CORRECTED_SUPPORT_ONLY="#F28E2B",NOT_SUPPORTED_AFTER_CORRECTION="#E15759")
par(mar=c(7,5,4,2));barplot(ms$direction_agree_both_proportion,names.arg=sub("PROGRAM_MODULE_","M",ms$program_module),col=cols[ms$corrected_module_status],ylim=c(0,1),ylab="Constituent proportion agreeing in both subtypes",main="Corrected support across 11 legacy modules");abline(h=.5,lty=2);dev.off()

metrics <- data.table(
  metric=c("corrected_gene_rows","eligible_celltype_contrasts","frozen20_n","frozen20_direction_both_n","frozen20_fdr_any_n","frozen20_fdr_both_n","corrected_shared_support_n","tier2_tested_pairs","tier2_corrected_shared_n","tier2_eope_enriched_n","tier2_lope_enriched_n","corrected_module_support_n","directional_module_n","not_supported_module_n"),
  value=c(nrow(genes),nrow(rbindlist(diagnostic_rows)),nrow(t1w),sum(t1w$direction_agrees_EOPE=="YES"&t1w$direction_agrees_LOPE=="YES"),sum(t1w$BH_FDR_frozen20_EOPE<.05|t1w$BH_FDR_frozen20_LOPE<.05),sum(t1w$BH_FDR_frozen20_EOPE<.05&t1w$BH_FDR_frozen20_LOPE<.05),sum(t1w$classification=="CORRECTED_SHARED_SUPPORT"),nrow(t2w),sum(t2w$classification=="CORRECTED_SHARED_PE"),sum(t2w$classification=="CORRECTED_EOPE_ENRICHED"),sum(t2w$classification=="CORRECTED_LOPE_ENRICHED"),sum(module_summary$corrected_module_status=="CORRECTED_ADMATI_SUPPORT"),sum(module_summary$corrected_module_status=="DIRECTIONAL_CORRECTED_SUPPORT_ONLY"),sum(module_summary$corrected_module_status=="NOT_SUPPORTED_AFTER_CORRECTION"))
)
write_csv(metrics,file.path(out,"corrected_analysis/phase2a2_metrics.csv"))
writeLines(capture.output(sessionInfo()),file.path(out,"phase2a2_session_info.txt"))
message("Phase 2A.2 corrected Route B analysis complete")
