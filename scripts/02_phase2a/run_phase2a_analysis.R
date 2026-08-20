# DEPRECATED / HISTORICAL ONLY: Phase 2A edgeR count-likelihood receiver inference is not authoritative. Use scripts/02_phase2a2/run_phase2a2_corrected_analysis.R.
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
setwd(root)
set.seed(20260809)

out <- "results/02_phase2a"
dir.create(file.path(out, "DE"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "programs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "regulons"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
interim <- "data/interim/phase2a"

write_csv <- function(x, path) fwrite(x, path, na = "", quote = TRUE)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

message("Reading frozen pseudobulk counts")
raw <- fread(file.path(interim, "admati_harmonized_pseudobulk_counts.csv"), check.names = FALSE)
genes <- raw[[1]]
mat <- as.matrix(raw[, -1])
storage.mode(mat) <- "integer"
counts <- rowsum(mat, group = genes, reorder = FALSE)
rm(raw, mat); invisible(gc())
message(sprintf("Collapsed %d rows to %d unique gene symbols", length(genes), nrow(counts)))

patients <- fread(file.path(out, "metadata/patient_registry.csv"))
elig <- fread(file.path(out, "metadata/pseudobulk_eligibility.csv"))
pbreg <- fread(file.path(out, "pseudobulk/pseudobulk_registry.csv"))

# Outcome-blind marker verification of the frozen lineage mapping.
marker_map <- list(
  EVT=c("HLA-G","MMP2","ITGA5"), VCT=c("EGFR","TP63","KRT7"), SCT=c("CGA","CGB3","ERVW-1"),
  ENDOTHELIAL=c("PECAM1","VWF","KDR"), PERICYTE=c("RGS5","CSPG4","PDGFRB"),
  VASCULAR_SMOOTH_MUSCLE=c("ACTA2","MYH11","TAGLN"), PLACENTAL_STROMAL=c("COL1A1","COL3A1","DCN"),
  HOFBAUER=c("C1QC","CD163","FOLR2"), MACROPHAGE=c("CD68","LST1","FCGR3A"),
  MONOCYTE=c("S100A8","S100A9","FCN1"), NK=c("NKG7","GNLY","KLRD1"), T_CELL=c("CD3D","CD3E","TRAC"),
  B_CELL=c("CD79A","MS4A1","CD37"), NEUTROPHIL=c("FCGR3B","CSF3R","S100A8"),
  IMMUNE_PROGENITOR_OR_PROLIFERATING=c("CD34","MKI67","STMN1")
)
celltypes <- names(marker_map)
pooled <- sapply(celltypes, function(ct) rowSums(counts[, grepl(paste0("__", ct, "$"), colnames(counts)), drop=FALSE]))
pooled_cpm <- t(t(pooled) / pmax(colSums(pooled), 1) * 1e6)
marker_qc <- rbindlist(lapply(celltypes, function(ct) {
  g <- intersect(marker_map[[ct]], rownames(pooled_cpm))
  target <- if (length(g)) mean(pooled_cpm[g, ct]) else 0
  other <- if (length(g)) median(rowMeans(pooled_cpm[g, setdiff(celltypes, ct), drop=FALSE])) else 0
  data.table(harmonized_annotation=ct, marker_genes_expected=length(marker_map[[ct]]), marker_genes_present=length(g), marker_expression_cpm=target, marker_enrichment_vs_other=target/(other+1e-8), marker_validation_status=ifelse(length(g)>=2 && target/(other+1e-8)>=1.5,"PASS","FLAG_REVIEW"))
}))
ann <- fread(file.path(out, "metadata/celltype_annotation_registry.csv"))
ann <- ann[, names(ann)[!grepl("^marker_(genes_expected|genes_present|expression_cpm|enrichment_vs_other|validation_status)(\\.|$)", names(ann))], with=FALSE]
ann <- merge(ann, marker_qc, by="harmonized_annotation", all.x=TRUE)
ann[, annotation_confidence := ifelse(marker_validation_status=="PASS", "PUBLISHED_HIGH; MARKER_CHECK_PASS", "PUBLISHED; MARKER_CHECK_FLAG_REVIEW")]
write_csv(ann, file.path(out, "metadata/celltype_annotation_registry.csv"))

read_gmt <- function(path, collection) {
  lines <- readLines(path, warn=FALSE)
  result <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed=TRUE)[[1]]
    unique(fields[-c(1,2)])
  })
  names(result) <- paste0(collection, "::", vapply(lines, function(line) strsplit(line, "\t", fixed=TRUE)[[1]][1], character(1)))
  result
}
resources <- "data/raw/phase2a_resources"
pathways <- c(
  read_gmt(file.path(resources,"h.all.v2026.1.Hs.symbols.gmt"),"HALLMARK"),
  read_gmt(file.path(resources,"c2.cp.reactome.v2026.1.Hs.symbols.gmt"),"REACTOME"),
  read_gmt(file.path(resources,"c5.go.bp.v2026.1.Hs.symbols.gmt"),"GOBP")
)
message(sprintf("Frozen gene-set universe: %d sets", length(pathways)))

net <- fread(file.path(resources, "collectri_human_genesymbols_omnipath.tsv"), select=c("source_genesymbol","target_genesymbol","consensus_stimulation","consensus_inhibition"))
net <- net[(consensus_stimulation == TRUE & consensus_inhibition == FALSE) | (consensus_stimulation == FALSE & consensus_inhibition == TRUE)]
net[, weight := ifelse(consensus_stimulation, 1, -1)]
net <- unique(net[, .(tf=source_genesymbol, target=target_genesymbol, weight)])

make_meta <- function(rows, contrast) {
  m <- merge(rows, patients[, .(patient_id, disease_status, pe_subtype_or_control_group, gestational_age_group, delivery_gestational_age_weeks, female_fetus, iugr)], by="patient_id", all.x=TRUE)
  if (contrast == "EOPE") m[, disease := factor(ifelse(pe_subtype_or_control_group=="EOPE","PE","CONTROL"), levels=c("CONTROL","PE"))]
  if (contrast == "LOPE") m[, disease := factor(ifelse(pe_subtype_or_control_group=="LOPE","PE","CONTROL"), levels=c("CONTROL","PE"))]
  if (contrast == "COMBINED_PE_SECONDARY") {
    m[, disease := factor(disease_status, levels=c("CONTROL","PE"))]
    m[, onset_stratum := factor(gestational_age_group, levels=c("EARLY","LATE"))]
  }
  m
}

influence_metrics <- function(pathway_genes, logcpm, meta, contrast, full_direction) {
  g <- intersect(pathway_genes, rownames(logcpm))
  if (length(g) < 5) return(list(effect=NA_real_, max_cook=NA_real_, max_patient=NA_character_, loo_prop=NA_real_, all_same=NA_character_))
  score <- colMeans(logcpm[g,,drop=FALSE])
  dat <- data.frame(score=score, disease=meta$disease, onset_stratum=if ("onset_stratum" %in% names(meta)) meta$onset_stratum else factor(rep("ONE",nrow(meta))), patient_id=meta$patient_id)
  form <- if (contrast=="COMBINED_PE_SECONDARY") score ~ onset_stratum + disease else score ~ disease
  fit <- lm(form, data=dat)
  effect <- unname(coef(fit)["diseasePE"])
  cooks <- cooks.distance(fit)
  effects <- vapply(seq_len(nrow(dat)), function(i) {
    f <- try(lm(form, data=dat[-i,,drop=FALSE]), silent=TRUE)
    if (inherits(f,"try-error")) return(NA_real_)
    unname(coef(f)["diseasePE"])
  }, numeric(1))
  valid <- is.finite(effects)
  same <- if (any(valid)) mean(sign(effects[valid]) == sign(full_direction)) else NA_real_
  list(effect=effect, max_cook=max(cooks,na.rm=TRUE), max_patient=dat$patient_id[which.max(cooks)], loo_prop=same, all_same=ifelse(isTRUE(all(sign(effects[valid])==sign(full_direction))),"YES","NO"))
}

run_regulons <- function(stat, celltype, contrast) {
  rbindlist(lapply(split(net, net$tf), function(x) {
    response <- unname(stat[x$target])
    ok <- is.finite(response)
    response <- response[ok]; weight <- x$weight[ok]
    n <- length(response)
    if (n < 10 || length(unique(weight)) < 2) return(NULL)
    xc <- weight-mean(weight); yc <- response-mean(response)
    sxx <- sum(xc^2); slope <- sum(xc*yc)/sxx
    intercept <- mean(response)-slope*mean(weight)
    rss <- sum((response-intercept-slope*weight)^2)
    se <- sqrt((rss/(n-2))/sxx); statistic <- slope/se
    data.table(celltype=celltype, contrast=contrast, tf=x$tf[1], target_n=n, activity_effect=slope, SE=se, statistic=statistic, P=2*pt(-abs(statistic),df=n-2))
  }), fill=TRUE)
}

contrast_order <- c("EOPE","LOPE","COMBINED_PE_SECONDARY")
de_summary <- list(); program_results <- list(); regulon_results <- list(); mds_rows <- list(); de_index <- list()
for (contrast_name in contrast_order) {
  contrast <- contrast_name
  eligible_ct <- sort(unique(elig[contrast==contrast_name & celltype_contrast_eligible=="YES", harmonized_annotation]))
  for (ct in eligible_ct) {
    rows <- elig[contrast==contrast_name & harmonized_annotation==ct & include_in_contrast=="YES"]
    meta <- make_meta(rows, contrast)
    cols <- rows$patient_id
    matrix_cols <- paste0(cols,"__",ct)
    x <- counts[, matrix_cols, drop=FALSE]
    colnames(x) <- cols
    keep <- rowSums(cpm(x) >= 1) >= 3
    y <- DGEList(x[keep,,drop=FALSE])
    y <- calcNormFactors(y, method="TMM")
    design <- if (contrast=="COMBINED_PE_SECONDARY") model.matrix(~ onset_stratum + disease, data=meta) else model.matrix(~ disease, data=meta)
    stopifnot(qr(design)$rank == ncol(design))
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmQLFit(y, design, robust=TRUE)
    coef_name <- grep("diseasePE", colnames(design), value=TRUE)
    test <- glmQLFTest(fit, coef=match(coef_name,colnames(design)))
    tab <- as.data.table(topTags(test,n=Inf,sort.by="none")$table, keep.rownames="gene")
    tab[, statistic := sign(logFC)*sqrt(pmax(F,0))]
    tab[, SE := fifelse(F>0,abs(logFC)/sqrt(pmax(F,0)),NA_real_)]
    tab[, `:=`(P=PValue, BH_FDR=FDR, n_PE=sum(meta$disease=="PE"), n_control=sum(meta$disease=="CONTROL"), contrast=contrast, celltype=ct, model_formula=ifelse(contrast=="COMBINED_PE_SECONDARY","~ onset_stratum + disease","~ disease"), SE_method="abs(logFC)/sqrt(QLF); descriptive 1-df equivalent")]
    setcolorder(tab,c("gene","logFC","SE","statistic","P","BH_FDR","logCPM","n_PE","n_control","contrast","celltype","model_formula","SE_method","F","PValue","FDR"))
    de_path <- file.path(out,"DE",paste0(contrast,"__",safe_name(ct),"_DE.csv"))
    write_csv(tab,de_path)
    de_index[[paste(contrast,ct)]] <- tab
    lcp <- cpm(y, log=TRUE, prior.count=2)
    mds <- plotMDS(lcp, plot=FALSE)
    d <- sqrt((mds$x-mean(mds$x))^2+(mds$y-mean(mds$y))^2)
    z <- (d-median(d))/(mad(d,constant=1)+1e-8)
    mds_rows[[length(mds_rows)+1]] <- data.table(contrast=contrast,celltype=ct,patient_id=meta$patient_id,disease=as.character(meta$disease),MDS1=mds$x,MDS2=mds$y,robust_distance_z=z,outlier_flag=ifelse(z>3.5,"FLAG_REVIEW","RETAIN"))
    de_summary[[length(de_summary)+1]] <- data.table(contrast=contrast,celltype=ct,n_PE=sum(meta$disease=="PE"),n_control=sum(meta$disease=="CONTROL"),eligible_gene_n=nrow(tab),DE_FDR05_n=sum(tab$BH_FDR<0.05),design_rank=qr(design)$rank,design_columns=ncol(design),full_rank="YES",primary_or_secondary=ifelse(contrast=="COMBINED_PE_SECONDARY","SECONDARY_SHARED_PROGRAM_SUPPORT","PRIMARY_STRATIFIED"),de_file=de_path)

    stats <- tab$statistic; names(stats) <- tab$gene; stats <- sort(stats[is.finite(stats)],decreasing=TRUE)
    pathway_index <- lapply(pathways, function(g) {
      ix <- match(g, names(stats), nomatch=0L)
      ix[ix>0L]
    })
    pathway_index <- pathway_index[lengths(pathway_index)>=15 & lengths(pathway_index)<=500]
    fg <- as.data.table(cameraPR(stats, pathway_index, inter.gene.cor=0.01, sort=FALSE), keep.rownames="pathway")
    fg[, size := NGenes]
    fg[, ES := vapply(pathway_index[pathway], function(ix) mean(stats[ix]), numeric(1))]
    fg[, NES := ES / sd(stats)]
    fg[Direction=="Down", NES := -abs(NES)]
    fg[Direction=="Up", NES := abs(NES)]
    fg[, P := PValue]
    fg[, BH_FDR := p.adjust(P,"BH")]
    fg[, c("collection","gene_set") := tstrsplit(pathway,"::",fixed=TRUE,keep=1:2)]
    fg[, direction := ifelse(NES>0,"UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM")]
    fg[, leadingEdge := "NOT_AVAILABLE_CAMERA_PR"]
    fg[, `:=`(celltype=ct,contrast=contrast,n_PE=sum(meta$disease=="PE"),n_control=sum(meta$disease=="CONTROL"),gene_set_release="MSigDB 2026.1.Hs",classification="UNCLASSIFIED_PENDING_SUBTYPE_MERGE",patient_score_effect=NA_real_,max_cooks_distance=NA_real_,max_influence_patient=NA_character_,loo_direction_same_proportion=NA_real_,all_loo_directions_same=NA_character_)]
    sigix <- which(fg$BH_FDR < 0.05)
    if (length(sigix)) for (i in sigix) {
      inf <- influence_metrics(pathways[[fg$pathway[i]]], lcp, meta, contrast, fg$NES[i])
      fg[i, `:=`(patient_score_effect=inf$effect,max_cooks_distance=inf$max_cook,max_influence_patient=inf$max_patient,loo_direction_same_proportion=inf$loo_prop,all_loo_directions_same=inf$all_same)]
    }
    program_results[[length(program_results)+1]] <- fg[, .(celltype,contrast,collection,gene_set,pathway,size,ES,NES,P,BH_FDR,direction,leadingEdge,n_PE,n_control,gene_set_release,classification,patient_score_effect,max_cooks_distance,max_influence_patient,loo_direction_same_proportion,all_loo_directions_same)]
    reg <- run_regulons(stats,ct,contrast)
    if (nrow(reg)) { reg[, BH_FDR := p.adjust(P,"BH")]; reg[, direction := ifelse(activity_effect>0,"UP_ACTIVITY_IN_PE","DOWN_ACTIVITY_IN_PE")]; reg[, `:=`(network="CollecTRI_OmniPath_2026-08-09",method="ULM_ON_SIGNED_DE_STATISTICS",exploratory="YES")]; regulon_results[[length(regulon_results)+1]] <- reg }
    message(sprintf("%s %s: %d genes, %d DE FDR<.05, %d programs FDR<.05",contrast,ct,nrow(tab),sum(tab$BH_FDR<.05),sum(fg$BH_FDR<.05)))
  }
}

de_summary_dt <- rbindlist(de_summary)
write_csv(de_summary_dt,file.path(out,"DE/celltype_DE_summary.csv"))
mds_dt <- rbindlist(mds_rows)
write_csv(mds_dt,file.path(out,"qc/pseudobulk_mds_coordinates.csv"))
programs <- rbindlist(program_results,fill=TRUE)

# Frozen classification based only on independently run EOPE and LOPE contrasts.
wide <- dcast(programs[contrast %in% c("EOPE","LOPE")], celltype+collection+gene_set+pathway ~ contrast, value.var=c("NES","BH_FDR","all_loo_directions_same"))
wide[, classification := "NOT_SIGNIFICANT"]
wide[is.finite(BH_FDR_EOPE) & is.finite(BH_FDR_LOPE) & BH_FDR_EOPE<0.05 & BH_FDR_LOPE<0.05 & sign(NES_EOPE)==sign(NES_LOPE), classification := "SHARED_PE"]
wide[is.finite(BH_FDR_EOPE) & BH_FDR_EOPE<0.05 & (!is.finite(BH_FDR_LOPE) | BH_FDR_LOPE>=0.10), classification := "EOPE_ENRICHED"]
wide[is.finite(BH_FDR_LOPE) & BH_FDR_LOPE<0.05 & (!is.finite(BH_FDR_EOPE) | BH_FDR_EOPE>=0.10), classification := "LOPE_ENRICHED"]
wide[((is.finite(BH_FDR_EOPE)&BH_FDR_EOPE<0.05)|(is.finite(BH_FDR_LOPE)&BH_FDR_LOPE<0.05)) & is.finite(NES_EOPE)&is.finite(NES_LOPE)&sign(NES_EOPE)!=sign(NES_LOPE), classification := "UNSTABLE"]
wide[((is.finite(BH_FDR_EOPE)&BH_FDR_EOPE<0.05 & is.finite(BH_FDR_LOPE)&BH_FDR_LOPE>=0.05&BH_FDR_LOPE<0.10) | (is.finite(BH_FDR_LOPE)&BH_FDR_LOPE<0.05 & is.finite(BH_FDR_EOPE)&BH_FDR_EOPE>=0.05&BH_FDR_EOPE<0.10)), classification := "UNSTABLE"]
classmap <- wide[,.(celltype,pathway,classification)]
programs <- merge(programs[,setdiff(names(programs),"classification"),with=FALSE],classmap,by=c("celltype","pathway"),all.x=TRUE)
programs[contrast=="COMBINED_PE_SECONDARY",classification:="SECONDARY_COMBINED_SUPPORT_ONLY"]
programs[is.na(classification),classification:="NOT_SIGNIFICANT"]
write_csv(programs,file.path(out,"programs/pe_cellstate_programs.csv"))
write_csv(wide[classification=="SHARED_PE"],file.path(out,"programs/shared_pe_programs.csv"))
write_csv(programs[classification=="EOPE_ENRICHED" & contrast=="EOPE"],file.path(out,"programs/eope_specific_programs.csv"))
write_csv(programs[classification=="LOPE_ENRICHED" & contrast=="LOPE"],file.path(out,"programs/lope_specific_programs.csv"))
regulons <- rbindlist(regulon_results,fill=TRUE)
write_csv(regulons,file.path(out,"regulons/cellstate_regulon_activity.csv"))

robust <- programs[contrast %in% c("EOPE","LOPE") & BH_FDR<0.05 & all_loo_directions_same=="YES"]
robust_counts <- unique(robust[,.(celltype,pathway,classification)])[,.(robust_program_n=.N),by=celltype]
shared_n <- unique(robust[classification=="SHARED_PE",.(celltype,pathway)])[, .N]
receiver_ge3 <- robust_counts[robust_program_n>=3,.N]
gate <- if (nrow(robust)==0) "NO_GO" else if (receiver_ge3>=2 && shared_n>=1) "GO_TO_PHASE2B_WITH_RESTRICTIONS" else "REMAIN_IN_PHASE2A"

qc <- rbindlist(list(
  data.table(metric="independent_patient_n",value=nrow(patients),detail="10 EOPE; 7 LOPE; 3 early controls; 6 late controls"),
  data.table(metric="cell_n",value=sum(patients$cell_count),detail="published processed UMI table"),
  data.table(metric="library_n",value=sum(patients$library_count),detail="31 libraries collapse to 26 donors"),
  data.table(metric="original_annotation_n",value=nrow(ann),detail="published labels"),
  data.table(metric="harmonized_annotation_n",value=uniqueN(ann$harmonized_annotation),detail="outcome-blind frozen mapping"),
  data.table(metric="marker_check_pass_n",value=marker_qc[marker_validation_status=="PASS",.N],detail="canonical marker pooled-count check"),
  data.table(metric="EOPE_eligible_celltype_n",value=de_summary_dt[contrast=="EOPE",.N],detail="10 cases versus 3 early controls maximum"),
  data.table(metric="LOPE_eligible_celltype_n",value=de_summary_dt[contrast=="LOPE",.N],detail="7 cases versus 6 late controls maximum"),
  data.table(metric="significant_program_rows",value=programs[contrast %in% c("EOPE","LOPE") & BH_FDR<0.05,.N],detail="subtype-stratified rows; paired SHARED program counted twice"),
  data.table(metric="robust_program_rows",value=nrow(robust),detail="all valid leave-one-patient-out directions retained"),
  data.table(metric="robust_shared_program_pairs",value=shared_n,detail="unique celltype x gene set"),
  data.table(metric="phase2a_gate",value=gate,detail="frozen gate; early-control n=3 and count-layer mismatch force restrictions at best")
),fill=TRUE)
write_csv(qc,file.path(out,"qc/phase2a_qc_summary.csv"))

risks <- data.table(
  risk_id=c("P2A-001","P2A-002","P2A-003","P2A-004","P2A-005","P2A-006"),
  severity=c("HIGH","HIGH","HIGH","MEDIUM","MEDIUM","HIGH"),
  status=c("OPEN","OPEN","MITIGATED","OPEN","MITIGATED","OPEN"),
  risk=c("Only three gestational-age-compatible early controls","Published total_molecules does not equal count-matrix column sums","Cells/libraries could be mistaken for biological replicates","IUGR is present only in PE and cannot be separated cleanly from disease","Published annotations could encode over-clustered states","Single primary cohort cannot establish cross-dataset reproducibility"),
  impact=c("EOPE estimates and leave-one-out stability are fragile","Count-layer provenance is incomplete; matrix column sums used for normalization","Inflated significance if ignored","EOPE programs may partly reflect FGR biology","Subtype cell states might be compositionally sensitive","Phase 2B bulk support is necessary before receiver programs are used downstream"),
  mitigation=c("Patient-level QL models; explicit LOPO; restricted gate","No imputation/rescaling; flag mismatch; use integer count table as published","donorID is sole inferential unit; multiple libraries collapsed","Report IUGR confounding; do not claim PE-specific causality","Collapse to 15 lineage-level annotations and validate markers without outcomes","Freeze programs now; test them directionally across six bulk cohorts only in Phase 2B"),
  source_url=c("results/02_phase2a/metadata/patient_registry.csv","https://doi.org/10.6084/m9.figshare.23264102.v1","results/02_phase2a/metadata/patient_registry.csv","results/02_phase2a/metadata/patient_registry.csv","https://doi.org/10.1016/j.medj.2023.07.005","docs/DATASET_AUDIT_PHASE0B_REPORT.md")
)
write_csv(risks,file.path(out,"qc/phase2a_risk_flags.csv"))

# Analytical QC figures only; no therapeutic interpretation.
png(file.path(out,"figures/A_cells_per_patient.png"),1200,700,res=130)
cols <- c(CONTROL="#4C78A8",PE="#E45756")[patients$disease_status]
barplot(patients$cell_count,names.arg=patients$patient_id,las=2,col=cols,main="Published cells per patient",ylab="Cells",cex.names=.65); legend("topright",legend=c("Control","PE"),fill=c("#4C78A8","#E45756"),bty="n"); dev.off()

heat <- dcast(pbreg,patient_id~harmonized_annotation,value.var="cell_count"); hm <- as.matrix(heat[,-1]); rownames(hm)<-heat$patient_id
png(file.path(out,"figures/B_cells_patient_celltype_heatmap.png"),1400,1000,res=130)
image(t(log10(hm+1)[nrow(hm):1,]),axes=FALSE,col=hcl.colors(40,"YlOrRd"),main="log10(cells + 1) per patient × cell type"); axis(1,at=seq(0,1,length.out=ncol(hm)),labels=colnames(hm),las=2,cex.axis=.65); axis(2,at=seq(0,1,length.out=nrow(hm)),labels=rev(rownames(hm)),las=2,cex.axis=.65); dev.off()

png(file.path(out,"figures/C_pseudobulk_library_sizes.png"),1100,700,res=130)
boxplot(log10(pseudobulk_library_umi+1)~harmonized_annotation,data=pbreg,las=2,col="#72B7B2",ylab="log10(count-matrix library UMI + 1)",main="Patient pseudobulk library sizes",cex.axis=.65); abline(h=log10(10000+1),lty=2,col="red"); dev.off()

top_ct <- de_summary_dt[contrast %in% c("EOPE","LOPE"),.(n=sum(n_PE+n_control)),by=celltype][order(-n)][1:min(6,.N),celltype]
plot_mds <- mds_dt[contrast %in% c("EOPE","LOPE") & celltype %in% top_ct]
png(file.path(out,"figures/D_patient_pseudobulk_MDS.png"),1400,900,res=130); par(mfrow=c(2,3),mar=c(4,4,3,1))
for (ct in top_ct) { z<-plot_mds[celltype==ct]; plot(z$MDS1,z$MDS2,col=ifelse(z$disease=="PE","#E45756","#4C78A8"),pch=19,xlab="MDS1",ylab="MDS2",main=ct); text(z$MDS1,z$MDS2,labels=z$patient_id,pos=3,cex=.45) }; dev.off()

pc <- unique(programs[contrast %in% c("EOPE","LOPE") & BH_FDR<0.05,.(celltype,pathway,classification)])[, .N,by=.(celltype,classification)]
png(file.path(out,"figures/E_significant_program_counts.png"),1200,750,res=130)
if (nrow(pc)) { widepc<-dcast(pc,celltype~classification,value.var="N",fill=0); m<-as.matrix(widepc[,-1]); rownames(m)<-widepc$celltype; barplot(t(m),beside=FALSE,las=2,col=hcl.colors(ncol(m),"Set2"),main="Significant subtype-stratified program classifications",ylab="Unique gene sets"); legend("topright",legend=colnames(m),fill=hcl.colors(ncol(m),"Set2"),cex=.7,bty="n") } else plot.new(); dev.off()

top_shared <- head(unique(programs[classification=="SHARED_PE" & contrast %in% c("EOPE","LOPE"),.(celltype,pathway,collection,gene_set)]),30)
png(file.path(out,"figures/F_shared_program_heatmap.png"),1400,900,res=130)
if (nrow(top_shared)) { zz<-merge(programs[contrast %in% c("EOPE","LOPE"),.(celltype,pathway,contrast,NES)],top_shared,by=c("celltype","pathway")); zz[,label:=paste(celltype,gene_set,sep=" | ")]; w<-dcast(zz,label~contrast,value.var="NES"); mm<-as.matrix(w[,-1]); rownames(mm)<-w$label; image(t(mm[nrow(mm):1,,drop=FALSE]),axes=FALSE,col=hcl.colors(51,"Blue-Red 3",rev=TRUE),main="Shared PE programs: subtype NES"); axis(1,at=seq(0,1,length.out=ncol(mm)),labels=colnames(mm),las=2); axis(2,at=seq(0,1,length.out=nrow(mm)),labels=rev(rownames(mm)),las=2,cex.axis=.45) } else { plot.new(); text(.5,.5,"No SHARED_PE programs under frozen criteria") }; dev.off()

writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,"qc/phase2a_session_info.txt"))
writeLines(gate,file.path(interim,"phase2a_gate.txt"))
message("PHASE2A_GATE=",gate)
