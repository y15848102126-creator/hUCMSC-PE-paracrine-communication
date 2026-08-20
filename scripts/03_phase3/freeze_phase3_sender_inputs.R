#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table);library(digest)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1],winslash="/",mustWork=TRUE) else normalizePath(".",winslash="/")
setwd(root)
out <- "results/03_phase3"
for(d in c("metadata","ligand_universe")) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
write_csv <- function(x,p) fwrite(x,p,na="",quote=TRUE)

cfg <- jsonlite::fromJSON("config/phase3_analysis.json",simplifyVector=FALSE)
stopifnot(cfg$outcome_inspection_status_at_freeze=="NO_HUCMSC_LIGAND_IDENTITIES_OR_EXPRESSION_OUTCOMES_INSPECTED")
for(p in names(cfg$upstream_sha256)) stopifnot(tolower(digest(file=p,algo="sha256",serialize=FALSE))==cfg$upstream_sha256[[p]])
manifest <- fread("data/raw/phase3/download_manifest.csv")
stopifnot(nrow(manifest)==48,all(file.exists(file.path("data/raw/phase3",manifest$filename))))
for(i in seq_len(nrow(manifest))) stopifnot(digest(file=file.path("data/raw/phase3",manifest$filename[i]),algo="sha256",serialize=FALSE)==manifest$sha256[i])

read_soft <- function(acc) readLines(gzfile(file.path("data/raw",paste0(acc,"_family.soft.gz")),encoding="UTF-8"),warn=FALSE)
sample_blocks <- function(acc){
  z<-read_soft(acc);idx<-grep("^\\^SAMPLE = ",z);ends<-c(idx[-1]-1,length(z));lapply(seq_along(idx),function(i)z[idx[i]:ends[i]])
}
field1 <- function(block,prefix){z<-sub(prefix,"",grep(paste0("^",prefix),block,value=TRUE));if(length(z))z[1]else NA_character_}
chars <- function(block,key){z<-grep("^!Sample_characteristics_ch1 = ",block,value=TRUE);z<-sub("^!Sample_characteristics_ch1 = ","",z);v<-z[startsWith(tolower(z),paste0(tolower(key),":"))];if(length(v))trimws(sub("^[^:]+:","",v[1]))else NA_character_}
barcode_n <- function(filename) length(readLines(gzfile(file.path("data/raw/phase3",filename)),warn=FALSE))

# Complete GSE182158 MSC context registry: only U01/U02 are sender-eligible UC donors.
b182 <- sample_blocks("GSE182158")
d182 <- rbindlist(lapply(b182,function(b){
  gsm<-sub("^\\^SAMPLE = ","",b[1]);id<-field1(b,"!Sample_title = ");tissue<-chars(b,"tissue")
  f<-manifest[grepl(paste0("^",gsm,"_.*barcodes.tsv.gz$"),filename),filename][1]
  data.table(dataset="GSE182158",sample_id=gsm,donor_id=id,independence_group=paste0("GSE182158_",id),is_hucmsc_sender_donor=ifelse(tissue=="Umbilical cord","YES","NO_CONTEXT_ONLY"),tissue_source=tissue,tissue_definition=ifelse(tissue=="Umbilical cord","umbilical-cord-derived MSC; Wharton's jelly not explicitly resolved",paste0(tissue,"-derived MSC context")),passage="P1_OR_P2_PROTOCOL_DONOR_SPECIFIC_UNRESOLVED",condition="UNTREATED",cytokine="NONE",duration_hours=NA_real_,sex=ifelse(id=="U01","female",ifelse(id=="U02","female",NA_character_)),age=ifelse(id=="U01","28",ifelse(id=="U02","37",NA_character_)),platform="10x Chromium 3' v2; NovaSeq 6000",raw_count_availability="TRUE_CELLRANGER_COUNTS",normalized_availability="PUBLIC_PROCESSED_IN_PUBLICATION",cells_total=barcode_n(f),cells_untreated=barcode_n(f),cells_stimulated=0L,batch_structure="one Cell Ranger matrix per donor; tissue/source and protocol are inseparable",source_url=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gsm,"|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182158"))
}))

# Four GSE199071 HUCMSC donors; HUVEC samples are explicitly outside sender scope.
b199 <- sample_blocks("GSE199071")
d199 <- rbindlist(lapply(b199,function(b){
  title<-field1(b,"!Sample_title = ");if(!startsWith(toupper(title),"HUCMSC"))return(NULL)
  gsm<-sub("^\\^SAMPLE = ","",b[1]);id<-sub("_.*$","",title);f<-manifest[grepl(paste0("^",gsm,"_.*barcodes.tsv.gz$"),filename),filename][1];n<-barcode_n(f)
  data.table(dataset="GSE199071",sample_id=gsm,donor_id=id,independence_group=paste0("GSE199071_",id),is_hucmsc_sender_donor="YES",tissue_source="Umbilical cord",tissue_definition="human umbilical-cord MSC; exact Wharton's jelly versus whole-cord definition unresolved",passage=paste0("P",chars(b,"passage")),condition="UNTREATED",cytokine="NONE",duration_hours=NA_real_,sex=tolower(chars(b,"sex")),age=NA_character_,platform="10x Chromium; NextSeq 550",raw_count_availability="TRUE_CELLRANGER_COUNTS",normalized_availability="CELLRANGER_OUTPUT",cells_total=n,cells_untreated=n,cells_stimulated=0L,batch_structure="donor1 P3/female; donors2-4 P6/male, so passage and sex are partially confounded with donor",source_url=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",gsm,"|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199071"))
}),fill=TRUE)

# GSE117837 cell-level sample metadata, summarized by donor and frozen strata.
b117 <- sample_blocks("GSE117837")
s117 <- rbindlist(lapply(b117,function(b){
  tr<-tolower(chars(b,"treatment"));condition<-ifelse(grepl("stim|ifn|tnf",tr),"IFNG_TNFA_STIMULATED","UNTREATED")
  data.table(sample_id=sub("^\\^SAMPLE = ","",b[1]),donor_id=chars(b,"donor id"),passage=paste0("P",chars(b,"passage")),condition=condition)
}))
stopifnot(nrow(s117)==361,all(s117$donor_id%in%c("Donor1","Donor2")))
strata<-s117[,.(cells=.N,sample_ids=paste(sample_id,collapse=";")),by=.(donor_id,passage,condition)]
wide<-dcast(strata,donor_id+passage~condition,value.var="cells",fill=0)
if(!"IFNG_TNFA_STIMULATED"%in%names(wide))wide[,IFNG_TNFA_STIMULATED:=0L]
if(!"UNTREATED"%in%names(wide))wide[,UNTREATED:=0L]
wide[,`:=`(dataset="GSE117837",stratum_id=paste(donor_id,passage,sep="_"),valid_within_stratum_contrast=ifelse(UNTREATED>0&IFNG_TNFA_STIMULATED>0,"YES","NO_UNPAIRED"),cytokine="IFN-gamma 10 ng/mL + TNF-alpha 10 ng/mL",duration_hours=12,biological_replication_unit="DONOR_NOT_CELL_NOT_PASSAGE",source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837|https://pmc.ncbi.nlm.nih.gov/articles/PMC6506509/")]
setcolorder(wide,c("dataset","stratum_id","donor_id","passage","UNTREATED","IFNG_TNFA_STIMULATED","valid_within_stratum_contrast","cytokine","duration_hours","biological_replication_unit","source_url"))
d117<-s117[,.(cells_total=.N,cells_untreated=sum(condition=="UNTREATED"),cells_stimulated=sum(condition=="IFNG_TNFA_STIMULATED"),passage=paste(sort(unique(passage)),collapse=";")),by=donor_id]
d117[,`:=`(dataset="GSE117837",sample_id="CELL_LEVEL_GSM_SERIES",independence_group=paste0("GSE117837_",donor_id),is_hucmsc_sender_donor="YES",tissue_source="Umbilical cord",tissue_definition="umbilical-cord-derived MSC; exact Wharton's jelly versus whole-cord definition unresolved",condition="UNTREATED_AND_IFNG_TNFA_STIMULATED",cytokine="IFN-gamma 10 ng/mL + TNF-alpha 10 ng/mL",duration_hours=12,sex=NA_character_,age=NA_character_,platform="Fluidigm C1; HiSeq 2500",raw_count_availability="SUBMITTED_RAW_READ_COUNT_TABLE",normalized_availability="RPKM_DESCRIBED_NOT_USED",batch_structure="donor, passage and treatment strata; multiple passages from Donor2 are not independent donors",source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837|https://pmc.ncbi.nlm.nih.gov/articles/PMC6506509/")]
setcolorder(d117,names(d182))
donors<-rbindlist(list(d182,d199,d117),fill=TRUE);donors[,biological_unit:="DONOR"];setorder(donors,dataset,donor_id)
write_csv(donors,file.path(out,"metadata/hucmsc_donor_registry.csv"))
write_csv(wide,file.path(out,"metadata/licensing_strata_registry.csv"))

roles<-data.table(dataset=c("GSE182158","GSE199071","GSE117837"),confirmed_role=c("HUCMSC_ATLAS_CONTEXT","INDEPENDENT_HUCMSC_BASELINE_REPLICATION","HUCMSC_INFLAMMATORY_LICENSING"),independent_hucmsc_donor_n=c(2,4,2),context_non_uc_msc_donor_n=c(9,0,0),formal_phase3_use=c("baseline core plus donor-aware source context; never sufficient alone","independent baseline core replication","within-donor/passage licensing; untreated P0 is non-contrast context"),expression_layer=c("true Cell Ranger raw counts","true Cell Ranger raw counts","submitted cell-level raw read counts"),platform=c("10x v2 / NovaSeq 6000","10x / NextSeq 550","Fluidigm C1 / HiSeq 2500"),source_definition=c("umbilical cord; WJ versus whole cord unresolved","umbilical cord; WJ versus whole cord unresolved","umbilical-cord-derived; WJ versus whole cord unresolved"),role_confirmed_from_metadata="YES",source_url=c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182158|https://pmc.ncbi.nlm.nih.gov/articles/PMC8715893/","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199071","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837|https://pmc.ncbi.nlm.nih.gov/articles/PMC6506509/"))
write_csv(roles,file.path(out,"metadata/hucmsc_dataset_role_registry.csv"))

# Frozen ligand universe: primary membership is NicheNet only; OmniPath is annotation-only.
lr<-readRDS("data/raw/phase3/lr_network_human_21122021.rds");stopifnot(all(c("from","to")%in%names(lr)))
ligands<-sort(unique(as.character(lr$from)));stopifnot(length(ligands)>1000,!anyNA(ligands))
op<-fread("data/raw/phase3/omnipath_intercell_20260810.tsv",na.strings=c("","NA"));stopifnot(all(c("genesymbol","category","parent","database","transmitter","secreted")%in%names(op)))
truth<-function(x)tolower(as.character(x))%in%c("true","1","yes")
opa<-op[genesymbol%in%ligands,.(omnipath_row_n=.N,omnipath_ligand_annotation=ifelse(any(grepl("ligand",tolower(paste(category,parent)))),"YES","NO"),omnipath_transmitter_annotation=ifelse(any(truth(transmitter)),"YES","NO"),omnipath_secreted_annotation=ifelse(any(truth(secreted)),"YES","NO"),omnipath_extracellular_annotation=ifelse(any(truth(secreted)|grepl("extracellular|secreted",tolower(paste(category,parent)))),"YES","NO"),omnipath_databases=paste(sort(unique(database)),collapse=";")),by=genesymbol]
universe<-data.table(gene=ligands);universe<-merge(universe,opa,by.x="gene",by.y="genesymbol",all.x=TRUE)
universe[,`:=`(nichenet_lr_edge_n=vapply(gene,function(g)sum(lr$from==g),integer(1)),ligand_annotation_source="NicheNet lr_network_human_21122021.rds unique from",database_version="NicheNet human LR network 2021-12-21; Zenodo 7074291",mapping_status=ifelse(is.na(omnipath_row_n),"NICHENET_ONLY_NOT_OMNIPATH_MAPPED","NICHENET_AND_OMNIPATH_MAPPED"),protein_secretome_annotation=ifelse(omnipath_secreted_annotation=="YES","OMNIPATH_SECRETED_ANNOTATION",ifelse(is.na(omnipath_row_n),"NOT_MAPPED_TO_OMNIPATH","NO_OMNIPATH_SECRETED_ANNOTATION")),nichenet_resource_sha256=digest(file="data/raw/phase3/lr_network_human_21122021.rds",algo="sha256",serialize=FALSE),omnipath_snapshot_sha256=digest(file="data/raw/phase3/omnipath_intercell_20260810.tsv",algo="sha256",serialize=FALSE),freeze_status="FROZEN_BEFORE_EXPRESSION_OUTCOMES",source_url="https://zenodo.org/records/7074291|https://omnipathdb.org/intercell")]
for(j in c("omnipath_ligand_annotation","omnipath_transmitter_annotation","omnipath_secreted_annotation","omnipath_extracellular_annotation"))set(universe,which(is.na(universe[[j]])),j,"NOT_MAPPED")
universe[is.na(omnipath_row_n),omnipath_row_n:=0L];universe[is.na(omnipath_databases),omnipath_databases:=""]
setcolorder(universe,c("gene","ligand_annotation_source","database_version","mapping_status","nichenet_lr_edge_n","omnipath_row_n","omnipath_ligand_annotation","omnipath_transmitter_annotation","omnipath_secreted_annotation","omnipath_extracellular_annotation","protein_secretome_annotation","omnipath_databases","nichenet_resource_sha256","omnipath_snapshot_sha256","freeze_status","source_url"))
write_csv(universe,file.path(out,"ligand_universe/frozen_ligand_universe.csv"))
cat(sprintf("Phase 3 sender freeze complete: %d registry rows, %d hUCMSC donors, %d ligands, %d licensing strata\n",nrow(donors),sum(donors$is_hucmsc_sender_donor=="YES"),nrow(universe),nrow(wide)))
