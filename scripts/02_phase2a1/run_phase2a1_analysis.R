#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

script_arg <- commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))][1]
root <- normalizePath(file.path(dirname(sub("^--file=","",script_arg)),"../.."),winslash="/",mustWork=TRUE)
setwd(root)
set.seed(20260809)

out <- "results/02_phase2a1"
interim <- "data/interim/phase2a1"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"figures"),recursive=TRUE,showWarnings=FALSE)
dir.create(interim,recursive=TRUE,showWarnings=FALSE)
write_csv <- function(x,path) fwrite(x,path,na="",quote=TRUE)
adm_source <- "https://doi.org/10.6084/m9.figshare.23264102.v1|https://pubmed.ncbi.nlm.nih.gov/37572658/"
yang_source <- "https://doi.org/10.3389/fimmu.2023.1142273|https://github.com/JustMoveOnnn/preeclampsia"

# Read-only freeze checks. The discovery files are never written in this script.
shared <- fread("results/02_phase2a/programs/shared_pe_programs.csv")
stopifnot(nrow(shared)==20L,uniqueN(shared[,.(celltype,pathway)])==20L,all(shared$classification=="SHARED_PE"))
phase2a_hashes <- tools::md5sum(c(
  "results/02_phase2a/programs/shared_pe_programs.csv",
  "results/02_phase2a/programs/eope_specific_programs.csv",
  "results/02_phase2a/programs/lope_specific_programs.csv",
  "results/02_phase2a/metadata/pseudobulk_eligibility.csv"
))
fwrite(data.table(file=names(phase2a_hashes),md5=unname(phase2a_hashes)),file.path(interim,"phase2a_frozen_hashes.csv"))

read_gmt <- function(path,collection) {
  lines <- readLines(path,warn=FALSE)
  x <- lapply(lines,function(z){f<-strsplit(z,"\t",fixed=TRUE)[[1]];unique(f[-c(1,2)])})
  names(x) <- paste0(collection,"::",vapply(lines,function(z)strsplit(z,"\t",fixed=TRUE)[[1]][1],character(1)))
  x
}
resources <- "data/raw/phase2a_resources"
pathways <- c(
  read_gmt(file.path(resources,"h.all.v2026.1.Hs.symbols.gmt"),"HALLMARK"),
  read_gmt(file.path(resources,"c2.cp.reactome.v2026.1.Hs.symbols.gmt"),"REACTOME"),
  read_gmt(file.path(resources,"c5.go.bp.v2026.1.Hs.symbols.gmt"),"GOBP")
)
stopifnot(all(shared$pathway %in% names(pathways)))

patients <- fread("results/02_phase2a/metadata/patient_registry.csv")
elig <- fread("results/02_phase2a/metadata/pseudobulk_eligibility.csv")
raw <- fread("data/interim/phase2a/admati_harmonized_pseudobulk_counts.csv",check.names=FALSE)
genes <- raw[[1]]
count_matrix <- rowsum(as.matrix(raw[,-1]),group=genes,reorder=FALSE)
storage.mode(count_matrix) <- "integer"
rm(raw); invisible(gc())

exact_perm_p <- function(values,labels) {
  ok <- is.finite(values) & !is.na(labels)
  values <- values[ok]; labels <- labels[ok]
  ncase <- sum(labels=="PE"); n <- length(values)
  if (ncase<1 || ncase>=n) return(NA_real_)
  observed <- mean(values[labels=="PE"])-mean(values[labels=="CONTROL"])
  combos <- combn(n,ncase)
  total <- sum(values)
  perm <- (colSums(matrix(values[combos],nrow=ncase))-0) / ncase -
    (total-colSums(matrix(values[combos],nrow=ncase))) / (n-ncase)
  mean(abs(perm) >= abs(observed)-1e-12)
}

score_sets <- function(counts,program_rows) {
  y <- DGEList(counts)
  y <- calcNormFactors(y,method="TMM")
  logcpm <- cpm(y,log=TRUE,prior.count=2)
  rankmat <- apply(logcpm,2,function(z) rank(z,ties.method="average")/(length(z)+1)-0.5)
  if (is.null(dim(rankmat))) rankmat <- matrix(rankmat,ncol=1,dimnames=list(rownames(logcpm),colnames(logcpm)))
  ans <- lapply(program_rows$pathway,function(p){
    g <- intersect(pathways[[p]],rownames(rankmat))
    if (length(g)<10) return(rep(NA_real_,ncol(rankmat)))
    colMeans(rankmat[g,,drop=FALSE])
  })
  ans <- do.call(rbind,ans); rownames(ans)<-program_rows$pathway; colnames(ans)<-colnames(counts)
  list(scores=ans,available=vapply(program_rows$pathway,function(p)length(intersect(pathways[[p]],rownames(rankmat))),integer(1)))
}

program_score_rows <- list(); score_cache <- list()
for (contrast in c("EOPE","LOPE")) {
  contrast_name <- contrast
  case_label <- contrast; control_label <- if (contrast=="EOPE") "EARLY_CONTROL" else "LATE_CONTROL"
  for (ct in unique(shared$celltype)) {
    pgr <- shared[celltype==ct]
    er <- elig[contrast==contrast_name & harmonized_annotation==ct & include_in_contrast=="YES"]
    ids <- er$patient_id; x <- count_matrix[,paste0(ids,"__",ct),drop=FALSE]; colnames(x)<-ids
    sc <- score_sets(x,pgr)
    meta <- patients[match(ids,patient_id)]
    labels <- ifelse(meta$pe_subtype_or_control_group==case_label,"PE","CONTROL")
    for (i in seq_len(nrow(pgr))) {
      vals <- sc$scores[i,]
      raw_eff <- mean(vals[labels=="PE"])-mean(vals[labels=="CONTROL"])
      original_stat <- if (contrast=="EOPE") pgr$NES_EOPE[i] else pgr$NES_LOPE[i]
      dsign <- sign(original_stat)
      program_score_rows[[length(program_score_rows)+1]] <- data.table(
        contrast=contrast,celltype=ct,collection=pgr$collection[i],gene_set=pgr$gene_set[i],pathway=pgr$pathway[i],
        discovery_direction=ifelse(dsign>0,"UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM"),available_gene_n=sc$available[i],
        n_PE=sum(labels=="PE"),n_control=sum(labels=="CONTROL"),patient_score_effect=raw_eff,
        direction_aligned_effect=raw_eff*dsign,direction_agrees=ifelse(raw_eff*dsign>0,"YES","NO"),
        exact_permutation_P=exact_perm_p(vals,labels),interpretation="METHOD_CONCORDANCE_SENSITIVITY; SAME_DISCOVERY_COHORT",
        source_url=adm_source
      )
      score_cache[[paste(contrast,ct,pgr$pathway[i],sep="|")]] <- data.table(patient_id=ids,label=labels,score=vals,iugr=meta$iugr,delivery_mode=meta$delivery_mode,induction=meta$induction)
    }
  }
}
patient_scores <- rbindlist(program_score_rows)
patient_scores[,BH_FDR:=p.adjust(exact_permutation_P,"BH",n=20),by=contrast]
setcolorder(patient_scores,c("contrast","celltype","collection","gene_set","pathway","discovery_direction","available_gene_n","n_PE","n_control","patient_score_effect","direction_aligned_effect","direction_agrees","exact_permutation_P","BH_FDR","interpretation","source_url"))
write_csv(patient_scores,file.path(out,"patient_program_score_sensitivity.csv"))

# Clinical confounding audit.
cramers_v <- function(tab) {
  n <- sum(tab); if (n==0) return(NA_real_)
  suppressWarnings(sqrt(as.numeric(chisq.test(tab,correct=FALSE)$statistic)/n))
}
hedges_g <- function(x,y) {
  nx<-length(x);ny<-length(y);df<-nx+ny-2
  if(nx<2||ny<2||df<=0)return(NA_real_)
  sp<-sqrt(((nx-1)*var(x)+(ny-1)*var(y))/df); if(!is.finite(sp)||sp==0)return(NA_real_)
  J<-1-3/(4*df-1); J*(mean(x)-mean(y))/sp
}
clinical <- list()
for (contrast in c("EOPE","LOPE")) {
  case_group<-contrast;control_group<-if(contrast=="EOPE")"EARLY_CONTROL" else "LATE_CONTROL"
  d<-patients[pe_subtype_or_control_group %in% c(case_group,control_group)]
  d[,group:=ifelse(pe_subtype_or_control_group==case_group,"PE","CONTROL")]
  cats<-list(delivery_mode="delivery_mode",induction="induction",IUGR_FGR="iugr",fetal_sex_female="female_fetus")
  for(nm in names(cats)){
    v<-cats[[nm]];tab<-table(d$group,d[[v]])
    zero_support<-any(colSums(tab)>0 & apply(tab,2,function(z)any(z==0)))
    mincell<-if(length(tab))min(tab)else 0
    est<-if(zero_support)"NON_ESTIMABLE" else if(mincell<2||min(table(d$group))<5)"WEAKLY_ESTIMABLE" else "ESTIMABLE"
    pe_values<-d[group=="PE"][[v]];control_values<-d[group=="CONTROL"][[v]]
    clinical[[length(clinical)+1]]<-data.table(contrast=contrast,variable=nm,variable_type="CATEGORICAL",n_PE=sum(d$group=="PE"),n_control=sum(d$group=="CONTROL"),
      PE_summary=paste(names(table(pe_values)),as.integer(table(pe_values)),collapse=";"),
      control_summary=paste(names(table(control_values)),as.integer(table(control_values)),collapse=";"),
      overlap_check=ifelse(zero_support,"POSITIVITY_FAILURE","LEVEL_SUPPORT_IN_BOTH_GROUPS"),test="FISHER_EXACT_TWO_SIDED",P_value=fisher.test(tab)$p.value,
      effect_size_name="CRAMERS_V",effect_size=cramers_v(tab),estimability=est,rationale=ifelse(zero_support,"At least one observed level has zero support in a disease group; no adjustment model is permitted.",ifelse(est=="WEAKLY_ESTIMABLE","Overlap exists but sparse cells/small n limit adjustment.","Observed levels have usable cross-group support.")),source_url=adm_source)
  }
  for (nm in c("delivery_gestational_age_weeks","maternal_age_years")) {
    x<-as.numeric(d[group=="PE"][[nm]]);y<-as.numeric(d[group=="CONTROL"][[nm]]);x<-x[is.finite(x)];y<-y[is.finite(y)]
    overlap<-max(min(x),min(y))<=min(max(x),max(y))
    est<-if(!overlap)"NON_ESTIMABLE"else if(min(length(x),length(y))<5)"WEAKLY_ESTIMABLE"else"ESTIMABLE"
    clinical[[length(clinical)+1]]<-data.table(contrast=contrast,variable=ifelse(nm=="delivery_gestational_age_weeks","gestational_age","maternal_age"),variable_type="CONTINUOUS",n_PE=length(x),n_control=length(y),
      PE_summary=sprintf("mean=%.3f;range=%.3f-%.3f",mean(x),min(x),max(x)),control_summary=sprintf("mean=%.3f;range=%.3f-%.3f",mean(y),min(y),max(y)),
      overlap_check=ifelse(overlap,sprintf("RANGE_OVERLAP=%.3f-%.3f",max(min(x),min(y)),min(max(x),max(y))),"NO_RANGE_OVERLAP"),test="WILCOXON_RANK_SUM_EXACT_WHERE_FEASIBLE",
      P_value=suppressWarnings(wilcox.test(x,y,exact=TRUE)$p.value),effect_size_name="HEDGES_G_PE_MINUS_CONTROL",effect_size=hedges_g(x,y),estimability=est,
      rationale=ifelse(est=="WEAKLY_ESTIMABLE","Ranges overlap, but a group has fewer than five patients.",ifelse(est=="ESTIMABLE","Ranges overlap with at least five patients per group.","No range overlap; no adjustment model is permitted.")),source_url=adm_source)
  }
}
clinical_dt<-rbindlist(clinical,fill=TRUE)
write_csv(clinical_dt,file.path(out,"clinical_confounding_audit.csv"))

# IUGR restriction and delivery-context sensitivities, using only the frozen patient scores.
iugr_rows<-list();delivery_rows<-list()
for(i in seq_len(nrow(patient_scores))){
  z<-patient_scores[i];cache<-score_cache[[paste(z$contrast,z$celltype,z$pathway,sep="|")]]
  restricted<-cache[label=="CONTROL" | (label=="PE" & iugr==0)]
  pe_n<-sum(restricted$label=="PE");co_n<-sum(restricted$label=="CONTROL")
  if(pe_n>=3&&co_n>=3){eff<-mean(restricted[label=="PE"]$score)-mean(restricted[label=="CONTROL"]$score);p<-exact_perm_p(restricted$score,restricted$label);status<-"ESTIMABLE"}else{eff<-NA_real_;p<-NA_real_;status<-"NON_ESTIMABLE_MINIMUM_DONOR_N"}
  sgn<-ifelse(z$discovery_direction=="UP_IN_PE_PROGRAM",1,-1)
  iugr_rows[[i]]<-data.table(contrast=z$contrast,celltype=z$celltype,collection=z$collection,gene_set=z$gene_set,pathway=z$pathway,discovery_direction=z$discovery_direction,
    non_IUGR_PE_n=pe_n,control_n=co_n,full_patient_score_effect=z$patient_score_effect,non_IUGR_patient_score_effect=eff,
    direction_aligned_effect=eff*sgn,direction_agrees=ifelse(is.finite(eff),ifelse(eff*sgn>0,"YES","NO"),"NOT_ESTIMABLE"),
    attenuation_ratio=ifelse(is.finite(eff)&&z$patient_score_effect!=0,eff/z$patient_score_effect,NA_real_),program_score_effect_change=eff-z$patient_score_effect,
    exact_permutation_P=p,estimability=status,interpretation="FROZEN_PROGRAM_RESTRICTION_SENSITIVITY_ONLY",source_url=adm_source)
  for(ctx in c("C_SECTION_ONLY","INDUCTION_ADJUSTMENT")){
    if(ctx=="C_SECTION_ONLY"){
      sub<-cache[delivery_mode=="C_SECTION"];pn<-sum(sub$label=="PE");cn<-sum(sub$label=="CONTROL")
      if(pn>=3&&cn>=3){deff<-mean(sub[label=="PE"]$score)-mean(sub[label=="CONTROL"]$score);dp<-exact_perm_p(sub$score,sub$label);dest<-"ESTIMABLE"}else{deff<-NA_real_;dp<-NA_real_;dest<-"NON_ESTIMABLE_POSITIVITY_OR_MINIMUM_N"}
    }else{
      pn<-sum(cache$label=="PE");cn<-sum(cache$label=="CONTROL");deff<-NA_real_;dp<-NA_real_;dest<-"NON_ESTIMABLE_POSITIVITY_FAILURE"
    }
    delivery_rows[[length(delivery_rows)+1]]<-data.table(contrast=z$contrast,context=ctx,celltype=z$celltype,collection=z$collection,gene_set=z$gene_set,pathway=z$pathway,discovery_direction=z$discovery_direction,
      PE_n=pn,control_n=cn,context_patient_score_effect=deff,direction_aligned_effect=deff*sgn,direction_agrees=ifelse(is.finite(deff),ifelse(deff*sgn>0,"YES","NO"),"NOT_ESTIMABLE"),
      attenuation_ratio=ifelse(is.finite(deff)&&z$patient_score_effect!=0,deff/z$patient_score_effect,NA_real_),exact_permutation_P=dp,estimability=dest,
      rationale=ifelse(ctx=="INDUCTION_ADJUSTMENT","Induction lacks cross-class positivity; no regression adjustment is fit.",ifelse(dest=="ESTIMABLE","Both groups retain at least three C-section pregnancies.","No adequate same-delivery support in both groups.")),source_url=adm_source)
  }
}
iugr_dt<-rbindlist(iugr_rows);iugr_dt[,BH_FDR:=p.adjust(exact_permutation_P,"BH",n=20),by=contrast]
delivery_dt<-rbindlist(delivery_rows);delivery_dt[context=="C_SECTION_ONLY",BH_FDR:=p.adjust(exact_permutation_P,"BH",n=20),by=contrast]
write_csv(iugr_dt,file.path(out,"shared_program_iugr_sensitivity.csv"))
write_csv(delivery_dt,file.path(out,"shared_program_delivery_sensitivity.csv"))

# Audit and rename the custom cameraPR display statistic without changing P values or classifications.
programs<-fread("results/02_phase2a/programs/pe_cellstate_programs.csv")
camera_rows<-merge(programs[contrast%in%c("EOPE","LOPE"),.(celltype,contrast,collection,gene_set,pathway,size,ES,NES,PValue=P,BH_FDR,direction)],shared[,.(celltype,pathway)],by=c("celltype","pathway"))
camera_rows[,`:=`(cameraPR_direction=ifelse(direction=="UP_IN_PE_PROGRAM","Up","Down"),original_display_name="NES",official_cameraPR_NES="NO",
  statistic_formula="SIGNED_ABS(mean(member signed edgeR QL statistics) / sd(all estimable gene statistics)); sign forced to cameraPR Direction",
  replacement_display_name="SIGNED_MEAN_STATISTIC_SD_SCALED",SIGNED_MEAN_STATISTIC_SD_SCALED=NES,
  significance_status="PValue and BH_FDR preserved exactly; no reclassification",source_url="https://bioconductor.org/packages/release/bioc/html/limma.html|scripts/02_phase2a/run_phase2a_analysis.R")]
setorder(camera_rows,celltype,pathway,contrast)
write_csv(camera_rows,file.path(out,"cameraPR_statistic_audit.csv"))

# Membership overlap and connected-component modules within cell type and direction.
nodes<-copy(shared)
nodes[,frozen_direction:=ifelse(NES_EOPE>0,"UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM")]
pair_rows<-list()
for(i in seq_len(nrow(nodes)))for(j in seq_len(nrow(nodes))){
  comparable<-nodes$celltype[i]==nodes$celltype[j]&&nodes$frozen_direction[i]==nodes$frozen_direction[j]
  a<-unique(pathways[[nodes$pathway[i]]]);b<-unique(pathways[[nodes$pathway[j]]]);inter<-length(intersect(a,b));uni<-length(union(a,b))
  pair_rows[[length(pair_rows)+1]]<-data.table(program_i=nodes$pathway[i],program_j=nodes$pathway[j],celltype_i=nodes$celltype[i],celltype_j=nodes$celltype[j],direction_i=nodes$frozen_direction[i],direction_j=nodes$frozen_direction[j],
    comparable_within_celltype_direction=ifelse(comparable,"YES","NO"),gene_n_i=length(a),gene_n_j=length(b),intersection_gene_n=ifelse(comparable,inter,NA_integer_),union_gene_n=ifelse(comparable,uni,NA_integer_),
    jaccard_similarity=ifelse(comparable,inter/uni,NA_real_),overlap_graph_edge=ifelse(comparable&&i!=j&&inter/uni>=0.25,"YES","NO"),source_url="MSigDB 2026.1.Hs; data/raw/phase2a_resources")
}
redundancy<-rbindlist(pair_rows)
write_csv(redundancy,file.path(out,"shared_program_redundancy_matrix.csv"))
adj<-split(seq_len(nrow(nodes)),paste(nodes$celltype,nodes$frozen_direction,sep="|"));module_id<-rep(NA_character_,nrow(nodes));counter<-0
for(ix in adj){
  unseen<-ix
  while(length(unseen)){
    counter<-counter+1;component<-unseen[1];front<-component
    while(length(front)){
      v<-front[1];front<-front[-1]
      neigh<-ix[vapply(ix,function(k){a<-pathways[[nodes$pathway[v]]];b<-pathways[[nodes$pathway[k]]];length(intersect(a,b))/length(union(a,b))>=0.25},logical(1))]
      new<-setdiff(neigh,component);component<-union(component,new);front<-union(front,new)
    }
    module_id[component]<-sprintf("PROGRAM_MODULE_%02d",counter);unseen<-setdiff(unseen,component)
  }
}
nodes[,program_module:=module_id]
module_summary<-nodes[,.(celltype=first(celltype),direction=first(frozen_direction),original_gene_set_n=.N,members=paste(gene_set,collapse=";"),pathways=paste(pathway,collapse=";")),by=program_module]
# Labels are assigned from the observed component membership, after graph construction.
module_summary[,module_label:=fcase(
  program_module=="PROGRAM_MODULE_01","HOFBAUER_INTERFERON_SIGNALING_CORE",
  program_module=="PROGRAM_MODULE_02","HOFBAUER_ANTIVIRAL_GENOME_RESTRICTION",
  program_module=="PROGRAM_MODULE_03","HOFBAUER_INTERFERON_RESPONSE_HALLMARKS",
  program_module=="PROGRAM_MODULE_04","HOFBAUER_IFN_STIMULATED_HOST_RESPONSE",
  program_module=="PROGRAM_MODULE_05","MACROPHAGE_INTERFERON_ALPHA_BETA_SIGNALING",
  program_module=="PROGRAM_MODULE_06","PLACENTAL_STROMAL_INTERFERON_ALPHA_RESPONSE",
  program_module=="PROGRAM_MODULE_07","PLACENTAL_STROMAL_INTERFERON_ALPHA_BETA_SIGNALING",
  program_module=="PROGRAM_MODULE_08","PLACENTAL_STROMAL_IFN_STIMULATED_HOST_RESPONSE",
  program_module=="PROGRAM_MODULE_09","SCT_MITOCHONDRIAL_RESPIRATION_CORE",
  program_module=="PROGRAM_MODULE_10","SCT_OXIDATIVE_PHOSPHORYLATION_HALLMARK",
  program_module=="PROGRAM_MODULE_11","SCT_REPLICATION_ORIGIN_TRANSITION",
  default=paste0(celltype,"_",program_module))]
nodes<-merge(nodes,module_summary[,.(program_module,module_label)],by="program_module",all.x=TRUE)
modules<-nodes[,.(record_type="ORIGINAL_GENE_SET",program_module,module_label,celltype,frozen_direction,collection,gene_set,pathway,original_BH_FDR_EOPE=BH_FDR_EOPE,original_BH_FDR_LOPE=BH_FDR_LOPE,source_url)]
modules<-rbind(modules,module_summary[,.(record_type="PROGRAM_MODULE",program_module,module_label,celltype,direction,collection="",gene_set=members,pathway=pathways,original_BH_FDR_EOPE=NA_real_,original_BH_FDR_LOPE=NA_real_,source_url="MSigDB 2026.1.Hs; frozen Phase 2A shared programs")],fill=TRUE)
write_csv(modules,file.path(out,"shared_program_modules.csv"))

# Independent Yang LOPE placenta replication: outcome-blind marker reconstruction.
marker_sets<-list(
  VCTp=c("MKI67","TOP2A","UBE2C","CENPF"),VCT=c("EGFR","TP63","KRT7","KRT8","KRT18"),EVT=c("HLA-G","MMP2","ITGA5","TAC3"),
  SCT=c("CGA","CGB3","CGB5","ERVW-1","PSG3","SDC1"),DC=c("FCER1A","CD1C","CLEC10A"),MC=c("C1QC","CD163","FOLR2","CD68","LST1"),
  T=c("CD3D","CD3E","TRAC"),NK=c("NKG7","GNLY","KLRD1"),B=c("MS4A1","CD79A","CD37"),VSM=c("ACTA2","TAGLN","MYH11","RGS5"),
  FB=c("COL1A1","COL3A1","DCN","LUM"),EB=c("HBB","HBA1","HBA2","GYPA")
)
yang_dir<-"data/interim/phase2a1/yang_uncompressed"
files<-list.files(yang_dir,pattern="^pla[1-6]_.*\\.txt$",full.names=TRUE)
yang_pb<-list();yang_qc<-list()
for(p in sort(files)){
  sid<-sub("_.*$","",basename(p));group<-ifelse(sid%in%c("pla1","pla2","pla3"),"CONTROL","LOPE")
  message("Yang matrix ",sid)
  d<-fread(p,check.names=FALSE);g<-d[[1]];m<-as.matrix(d[,-1]);storage.mode(m)<-"integer";rm(d);invisible(gc())
  lib<-colSums(m);feat<-colSums(m>0);mt<-colSums(m[grepl("^MT-",g),,drop=FALSE]);qc<-lib>2000&lib<150000&(mt/pmax(lib,1))<0.25
  norm<-log1p(t(t(m[,qc,drop=FALSE])/pmax(lib[qc],1))*1e4)
  scores<-sapply(marker_sets,function(ms){ix<-match(ms,g,nomatch=0);ix<-ix[ix>0];if(length(ix))colMeans(norm[ix,,drop=FALSE])else rep(NA_real_,ncol(norm))})
  detected<-sapply(marker_sets,function(ms){ix<-match(ms,g,nomatch=0);ix<-ix[ix>0];if(length(ix))colSums(m[ix,qc,drop=FALSE]>0)else rep(0,ncol(norm))})
  best<-max.col(scores,ties.method="first");sorted<-t(apply(scores,1,sort,decreasing=TRUE));margin<-sorted[,1]-sorted[,2];labs<-colnames(scores)[best]
  pass<-detected[cbind(seq_len(nrow(detected)),best)]>=2&margin>=0.05;labs[!pass]<-"UNRESOLVED"
  yang_qc[[sid]]<-data.table(sample_id=sid,group=group,input_cell_n=ncol(m),qc_cell_n=sum(qc),resolved_cell_n=sum(labs!="UNRESOLVED"),SCT_cell_n=sum(labs=="SCT"),FB_cell_n=sum(labs=="FB"),MC_cell_n=sum(labs=="MC"))
  for(label in c("SCT","FB")){
    ix<-which(labs==label)
    if(length(ix)){
      v<-rowSums(m[,which(qc)[ix],drop=FALSE]);names(v)<-g
      yang_pb[[paste(sid,label,sep="|")]]<-v
    }
  }
  rm(m,norm,scores,detected);invisible(gc())
}
allg<-sort(unique(unlist(lapply(yang_pb,names),use.names=FALSE)))
yang_counts<-matrix(0L,nrow=length(allg),ncol=length(yang_pb),dimnames=list(allg,names(yang_pb)))
for(nm in names(yang_pb))yang_counts[names(yang_pb[[nm]]),nm]<-yang_pb[[nm]]
fwrite(rbindlist(yang_qc),file.path(interim,"yang_annotation_qc.csv"))
yang_results<-list()
for(i in seq_len(nrow(shared))){
  z<-shared[i];published_label<-if(z$celltype=="SCT")"SCT"else if(z$celltype=="PLACENTAL_STROMAL")"FB"else if(z$celltype%in%c("HOFBAUER","MACROPHAGE"))"MC"else"NOT_AVAILABLE"
  map_status<-if(z$celltype=="SCT")"DIRECT_MATCH"else if(z$celltype=="PLACENTAL_STROMAL")"APPROXIMATE_MATCH"else"NOT_MAPPABLE"
  if(map_status=="NOT_MAPPABLE"){
    yang_results[[i]]<-data.table(celltype=z$celltype,yang_published_annotation=published_label,mapping_status=map_status,collection=z$collection,gene_set=z$gene_set,pathway=z$pathway,admati_direction=ifelse(z$NES_LOPE>0,"UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM"),available_gene_n=NA_integer_,LOPE_n=3L,control_n=3L,minimum_cells_per_donor=NA_integer_,patient_score_effect=NA_real_,direction_aligned_effect=NA_real_,direction_agrees="NOT_EVALUABLE",exact_permutation_P=NA_real_,BH_FDR=NA_real_,replication_interpretation="Yang MC cannot distinguish Hofbauer from other macrophage/DC states; no forced equivalence.",source_url=yang_source)
    next
  }
  cols<-grep(paste0("\\|",published_label,"$"),colnames(yang_counts),value=TRUE);ids<-sub("\\|.*$","",cols);counts<-yang_counts[,cols,drop=FALSE];colnames(counts)<-ids
  cq<-rbindlist(yang_qc);cellcol<-paste0(published_label,"_cell_n");mincells<-min(cq[[cellcol]])
  if(length(cols)!=6||mincells<20){eff<-pval<-NA_real_;agree<-"NOT_EVALUABLE_MINIMUM_CELL_OR_DONOR_N";avail<-length(intersect(pathways[[z$pathway]],rownames(counts)))}else{
    sc<-score_sets(counts,z);vals<-as.numeric(sc$scores[1,]);lab<-ifelse(ids%in%c("pla4","pla5","pla6"),"PE","CONTROL");eff<-mean(vals[lab=="PE"])-mean(vals[lab=="CONTROL"]);pval<-exact_perm_p(vals,lab);agree<-ifelse(eff*sign(z$NES_LOPE)>0,"YES","NO");avail<-sc$available[1]
  }
  yang_results[[i]]<-data.table(celltype=z$celltype,yang_published_annotation=published_label,mapping_status=map_status,collection=z$collection,gene_set=z$gene_set,pathway=z$pathway,admati_direction=ifelse(z$NES_LOPE>0,"UP_IN_PE_PROGRAM","DOWN_IN_PE_PROGRAM"),available_gene_n=avail,LOPE_n=3L,control_n=3L,minimum_cells_per_donor=mincells,patient_score_effect=eff,direction_aligned_effect=eff*sign(z$NES_LOPE),direction_agrees=agree,exact_permutation_P=pval,BH_FDR=NA_real_,replication_interpretation=ifelse(map_status=="DIRECT_MATCH","Direct SCT ontology match; same frozen rank-score method.","Placental FB is an approximate stromal match; interpret separately."),source_url=yang_source)
}
yang_dt<-rbindlist(yang_results,fill=TRUE);yang_dt[,BH_FDR:=p.adjust(exact_permutation_P,"BH",n=20)]
write_csv(yang_dt,file.path(out,"yang_lope_replication.csv"))

# Admati paper reconciliation metrics and anchor findings.
de_summary<-fread("results/02_phase2a/DE/celltype_DE_summary.csv")
recon<-de_summary[contrast%in%c("EOPE","LOPE"),.(celltype_n=.N,eligible_gene_tests=sum(eligible_gene_n),DE_FDR05_n=sum(DE_FDR05_n)),by=contrast]
for(con in c("EOPE","LOPE")){
  fl<-list.files("results/02_phase2a/DE",pattern=paste0("^",con,"__.*_DE.csv$"),full.names=TRUE)
  dd<-rbindlist(lapply(fl,fread),fill=TRUE)
  recon[contrast==con,`:=`(nominal_P05_n=sum(dd$P<0.05),median_abs_logFC=median(abs(dd$logFC)),q90_abs_logFC=quantile(abs(dd$logFC),0.9))]
}
anchors<-rbindlist(lapply(c("EOPE","LOPE"),function(con){
  rbindlist(lapply(c("SCT","HOFBAUER","MACROPHAGE","PLACENTAL_STROMAL"),function(ct){
    p<-file.path("results/02_phase2a/DE",paste0(con,"__",ct,"_DE.csv"));if(!file.exists(p))return(NULL);d<-fread(p);d[gene%in%c("FLT1","PGF","ISG15","IFIT1","MX1","STAT1"),.(contrast=con,celltype=ct,gene,logFC,P,BH_FDR)]
  }),fill=TRUE)
}),fill=TRUE)
fwrite(recon,file.path(interim,"admati_early_late_reconciliation.csv"));fwrite(anchors,file.path(interim,"admati_anchor_findings.csv"))

# Count-layer provenance: public analysis scripts normalize each cell to 10,000 and apply ceiling.
count_audit<-fread(file.path(interim,"admati_count_layer_cell_audit.csv"))
loader<-readLines("data/raw/phase2a1/admati_2023/admati_load_sc_PE_data_and_save_v1.m",warn=FALSE)
cap_line<-grep("tot_mol\\(tot_mol>5e4\\) = 5e4",loader,value=TRUE)
normalization_scripts<-list.files("data/raw/phase2a1/admati_2023",pattern="\\.m$",full.names=TRUE)
norm_lines<-unlist(lapply(normalization_scripts,function(p)grep("data = ceil\\(data\\./repmat\\(sum\\(data\\),length\\(data\\(:,1\\)\\),1\\)\\*10e3\\)",readLines(p,warn=FALSE),value=TRUE)))
count_resolved<-length(norm_lines)>0&&all(count_audit$consistent_with_ceil_to_10000=="YES")
count_summary<-data.table(cell_n=nrow(count_audit),published_sum=sum(count_audit$published_total_molecules),direct_matrix_sum=sum(count_audit$direct_matrix_sum),direct_sum_min=min(count_audit$direct_matrix_sum),direct_sum_max=max(count_audit$direct_matrix_sum),raw_total_exact_n=sum(count_audit$exact_raw_total=="YES"),ceiling_10000_consistent_n=sum(count_audit$consistent_with_ceil_to_10000=="YES"),public_loader_cap_line=cap_line,public_normalization_code=paste(unique(trimws(norm_lines)),collapse=" | "),resolution=ifelse(count_resolved,"RESOLVED_EXPRESSION_IS_CEILED_LIBRARY_SIZE_NORMALIZED_NOT_RAW_UMI","UNRESOLVED"),source_url="https://github.com/zeiselamit/PE_2023")
fwrite(count_summary,file.path(interim,"count_layer_provenance_summary.csv"))

# Risk register and compact analytical figures.
risks<-data.table(
  risk_id=paste0("P2A1-",sprintf("%03d",1:8)),severity=c("HIGH","HIGH","HIGH","MEDIUM","HIGH","MEDIUM","MEDIUM","MEDIUM"),
  status=c("OPEN","OPEN","OPEN",ifelse(count_resolved,"RESOLVED_DEFINITION_CRITICAL_METHOD_RISK","OPEN"),"OPEN","MITIGATED","OPEN","OPEN"),
  risk=c("EOPE delivery mode lacks C-section controls","Induction lacks disease-crossed positivity","IUGR occurs only among PE cases","total_molecules differs from uncapped matrix sums","Yang annotation objects/cell barcode labels are not public","cameraPR-derived NES name is misleading","Yang replication has only three pregnancies per group","Gene-set redundancy inflates the apparent count of biological themes"),
  impact=c("EOPE programs cannot be separated from delivery context in this cohort","Induction-adjusted disease effects are not identifiable","Full disease effects may include FGR biology","The Phase 2A edgeR count likelihood was applied to ceiled 10k-normalized values rather than raw UMI counts","Published 12-class annotation cannot be exactly reattached; marker-rule reconstruction is necessary","Readers could mistake a custom standardized mean statistic for GSEA NES","Only large directional effects can be detected","Twenty gene sets represent fewer nonredundant modules"),
  mitigation=c("Classify NON_ESTIMABLE; do not fit adjustment","No induction regression or forced subset","Non-IUGR restriction tests only frozen programs","Preserve Phase 2A unchanged; do not interpret its count-model FDR as definitive; revise receiver inference using a method appropriate for normalized expression or recover raw counts","Outcome-blind marker rules; DIRECT/APPROXIMATE/NOT_MAPPABLE statuses; no Hofbauer forcing","Rename to SIGNED_MEAN_STATISTIC_SD_SCALED while preserving P/FDR","Exact patient-level inference and directional theme interpretation only","Retain all 20 originals and add Jaccard graph modules"),
  source_url=c(adm_source,adm_source,adm_source,"https://github.com/zeiselamit/PE_2023/blob/main/load_sc_PE_data_and_save_v1.m",yang_source,"scripts/02_phase2a/run_phase2a_analysis.R",yang_source,"MSigDB 2026.1.Hs")
)
write_csv(risks,file.path(out,"phase2a1_risk_flags.csv"))

png(file.path(out,"figures","A_clinical_overlap.png"),1200,750,res=130)
par(mfrow=c(1,2),mar=c(8,4,3,1));for(con in c("EOPE","LOPE")){d<-patients[pe_subtype_or_control_group%in%c(con,ifelse(con=="EOPE","EARLY_CONTROL","LATE_CONTROL"))];d[,group:=ifelse(pe_subtype_or_control_group==con,"PE","CONTROL")];barplot(t(table(d$group,d$delivery_mode)),beside=TRUE,col=c("#4C78A8","#E45756"),las=2,main=paste(con,"delivery mode"),ylab="Pregnancies");legend("topright",legend=rownames(t(table(d$group,d$delivery_mode))),fill=c("#4C78A8","#E45756"),bty="n")};dev.off()
png(file.path(out,"figures","B_shared_program_stress_test.png"),1200,750,res=130)
tab<-rbind(patient_scores[,.(test="Full patient score",credible=sum(direction_agrees=="YES"),not_credible=sum(direction_agrees!="YES")),by=contrast],iugr_dt[,.(test="Non-IUGR",credible=sum(direction_agrees=="YES"),not_credible=sum(direction_agrees!="YES")),by=contrast]);m<-as.matrix(dcast(tab,test~contrast,value.var="credible")[,-1]);rownames(m)<-unique(tab$test);barplot(t(m),beside=TRUE,col=c("#59A14F","#F28E2B"),main="Frozen shared programs retaining discovery direction",ylab="Programs",legend.text=colnames(m));dev.off()
png(file.path(out,"figures","C_yang_replication.png"),1200,750,res=130)
yr<-yang_dt[is.finite(direction_aligned_effect)];cols<-ifelse(yr$direction_agrees=="YES","#59A14F","#E15759");short_names<-c("IFN-alpha\nhallmark","IFN-alpha/beta\nsignaling","IFN-stimulated\nhost response");par(mar=c(7,5,4,2));barplot(yr$direction_aligned_effect,names.arg=short_names,las=1,col=cols,main="Yang LOPE placenta: approximate stromal replication",ylab="Direction-aligned patient-score effect",cex.names=.9);abline(h=0,lty=2);legend("topright",c("Agrees","Opposes"),fill=c("#59A14F","#E15759"),bty="n");dev.off()
png(file.path(out,"figures","D_program_module_overlap.png"),1200,750,res=130)
ms<-module_summary[order(program_module)];par(mar=c(6,5,4,2));barplot(ms$original_gene_set_n,names.arg=sub("PROGRAM_MODULE_","M",ms$program_module),las=1,col=hcl.colors(nrow(ms),"Set2"),main="Eleven Jaccard-graph receiver modules",ylab="Frozen original gene sets",cex.names=.9);mtext("Module labels and memberships are reported in shared_program_modules.csv",side=1,line=4,cex=.75);dev.off()

metrics<-data.table(
  metric=c("shared_program_n","patient_score_direction_agree_EOPE","patient_score_direction_agree_LOPE","nonIUGR_direction_agree_EOPE","nonIUGR_direction_agree_LOPE","LOPE_Csection_direction_agree","yang_evaluable_program_n","yang_direction_agree_n","program_module_n","count_layer_resolution","EOPE_gene_DE_FDR05","LOPE_gene_DE_FDR05"),
  value=c(20,patient_scores[contrast=="EOPE",sum(direction_agrees=="YES")],patient_scores[contrast=="LOPE",sum(direction_agrees=="YES")],iugr_dt[contrast=="EOPE",sum(direction_agrees=="YES")],iugr_dt[contrast=="LOPE",sum(direction_agrees=="YES")],delivery_dt[contrast=="LOPE"&context=="C_SECTION_ONLY",sum(direction_agrees=="YES")],yang_dt[is.finite(patient_score_effect),.N],yang_dt[direction_agrees=="YES",.N],uniqueN(module_id),count_summary$resolution,recon[contrast=="EOPE",DE_FDR05_n],recon[contrast=="LOPE",DE_FDR05_n])
)
fwrite(metrics,file.path(interim,"phase2a1_metrics.csv"))
writeLines(sub("[[:space:]]+$","",capture.output(sessionInfo())),file.path(out,"phase2a1_session_info.txt"))
message("Phase 2A.1 analysis complete")
