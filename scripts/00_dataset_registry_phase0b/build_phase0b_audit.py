#!/usr/bin/env python3
"""Build Phase 0B feasibility registries from audited, cached primary evidence."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "phase0b"
OUT = ROOT / "results" / "00_dataset_audit_phase0b"
AUDIT_DATE = "2026-08-09"
CUTOFF = "2026-08-08"


def write_csv(name: str, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty registry: {name}")
    OUT.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    for row in rows:
        missing = set(fields) - set(row)
        extra = set(row) - set(fields)
        if missing or extra:
            raise ValueError(f"schema mismatch in {name}: missing={missing}, extra={extra}")
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def pe_row(**values: object) -> dict[str, object]:
    defaults: dict[str, object] = {
        "study_id": "", "title": "", "year": "", "publication_status": "",
        "pmid": "", "doi": "", "pe_subtype": "UNRESOLVED", "pe_donor_n": "UNRESOLVED",
        "control_donor_n": "UNRESOLVED", "other_donor_n": "0", "tissue": "",
        "compartment": "", "ga_case": "UNRESOLVED", "ga_control": "UNRESOLVED",
        "ga_matching": "UNRESOLVED", "platform": "UNRESOLVED", "repository_accession": "NOT_FOUND",
        "processed_data_availability": "NOT_FOUND", "raw_data_availability": "NOT_FOUND",
        "donor_identity_availability": "NO", "patient_level_pseudobulk_feasible": "NO",
        "possible_overlap": "NONE_IDENTIFIED", "recommended_role": "EXCLUDE",
        "exclusion_reason": "", "evidence_summary": "", "source_url": "",
        "source_accession": "", "source_pmid": "", "source_doi": "",
        "evidence_status": "AUDITED", "audit_date": AUDIT_DATE, "evidence_cutoff_date": CUTOFF,
    }
    defaults.update(values)
    return defaults


PE_ROWS = [
    pe_row(
        study_id="Admati_2023_FIGSHARE", title="Two distinct molecular faces of preeclampsia revealed by single-cell transcriptomics",
        year=2023, publication_status="PEER_REVIEWED", pmid="37572658", doi="10.1016/j.medj.2023.07.005",
        pe_subtype="EOPE n=10; LOPE n=7", pe_donor_n=17, control_donor_n=9,
        tissue="placenta", compartment="placental cotyledon", ga_case="approximately GW24-40; donor-level delivery_week is public",
        ga_control="3 preterm controls plus 6 term controls; donor-level delivery_week is public",
        ga_matching="SUBTYPE-AWARE CONTROLS AVAILABLE; model GA explicitly", platform="10x Genomics scRNA-seq",
        repository_accession="Figshare 23264102.v1; 23264165.v1", processed_data_availability="PUBLIC_PROCESSED: 86,752-cell expression matrix with metadata and annotations",
        raw_data_availability="NOT_FOUND", donor_identity_availability="YES: 26 stable donorID values and 31 library/sample labels",
        patient_level_pseudobulk_feasible="YES", possible_overlap="No external overlap identified; multiple libraries collapse to 26 donorID values",
        recommended_role="PRIMARY_PE_SCRNA", evidence_summary="Direct inspection of the 256.5 MB ZIP found 86,752 cell columns, 26 donor IDs (10 EOPE, 7 LOPE, 3 preterm controls, 6 term controls), cell annotations, delivery week, donor age, IUGR, fetal sex and delivery/treatment fields.",
        source_url="https://doi.org/10.1016/j.medj.2023.07.005|https://doi.org/10.6084/m9.figshare.23264102.v1|https://github.com/zeiselamit/PE_2023",
        source_accession="Figshare:23264102.v1;GitHub:zeiselamit/PE_2023", source_pmid="37572658", source_doi="10.1016/j.medj.2023.07.005",
    ),
    pe_row(
        study_id="PMID41472684_phs001886v6", title="Single-cell mapping of maternal-fetal cross-talk in preeclampsia",
        year=2025, publication_status="PREPRINT_AS_OF_CUTOFF", pmid="41472684", doi="10.21203/rs.3.rs-8254581/v1",
        pe_subtype="EOPE n=10; LOPE n=29", pe_donor_n=39, control_donor_n=39,
        tissue="maternal-fetal interface", compartment="decidua basalis and placental villi",
        ga_case="EOPE median 30.1 weeks; LOPE median 36.9 weeks (aggregate only)",
        ga_control="early-control median 30.5 weeks; late-control median 36.9 weeks (aggregate only)",
        ga_matching="AGGREGATE MATCHING CONFIRMED; individual mapping unavailable", platform="10x Genomics scRNA-seq plus maternal/fetal genotyping",
        repository_accession="phs001886.v6.p1; exact GEO accession NOT_FOUND",
        processed_data_availability="NOT_FOUND: no public matrix or exact GEO record located",
        raw_data_availability="CONTROLLED_RAW EXPECTED/CLAIMED, BUT v6 new sequence load not demonstrably available",
        donor_identity_availability="NO PUBLIC 78-patient case/control mapping", patient_level_pseudobulk_feasible="NO WITH PUBLIC FILES",
        possible_overlap="phs001886.v1-v5 are cumulative predecessors; v6 contains all prior public aliases",
        recommended_role="PENDING", evidence_summary="The preprint reports 39+39 singleton pregnancies, but public dbGaP v6 SSTR has 374 consented subjects/395 distinct samples and mixes maternal/fetal genotype records. New v6 records have no public BioSample/sequence mapping, and exact GEO title/author searches returned zero records. The 78 subjects cannot be reconstructed independently from public sample-level metadata.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12747280/|https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001886.v6.p1|https://doi.org/10.21203/rs.3.rs-8254581/v1",
        source_accession="phs001886.v6.p1;GEO=NOT_FOUND", source_pmid="41472684", source_doi="10.21203/rs.3.rs-8254581/v1",
    ),
    pe_row(
        study_id="GSE290578", title="Divide-and-conquer analysis reveals hidden immune cell influencers across the placenta in preeclampsia",
        year=2026, publication_status="PEER_REVIEWED", doi="10.1186/s13040-026-00556-y",
        pe_subtype="severe/early-onset PE", pe_donor_n=4, control_donor_n=4,
        tissue="placenta", compartment="paired maternal-fetal interface and deep placenta per pregnancy",
        ga_case="28+5 to 34+1 weeks", ga_control="38+2 to 39+6 weeks", ga_matching="FATAL FOR NAIVE DISEASE DISCOVERY: PE-preterm vs normal-term",
        platform="10x Genomics scRNA-seq; Illumina", repository_accession="GSE290578; PRJNA1228650",
        processed_data_availability="PUBLIC_PROCESSED: per-GSM Cell Ranger matrix triplets", raw_data_availability="PUBLIC_RAW: SRA",
        donor_identity_availability="YES: Norm1/4/6/7 and PePT2/3/5/8", patient_level_pseudobulk_feasible="YES, after re-annotation and paired-layer handling",
        possible_overlap="Analytically reuses phs001886.v1.p1 as negative validation; no subject overlap demonstrated",
        recommended_role="SCRNA_REPLICATION", evidence_summary="16 GSM are two paired tissue layers from eight women, not 16 biological replicates. Main analysis uses 3 PE+3 control after Norm1 QC failure; PePT8 is an independent within-study validation donor. Public matrices lack a stand-alone cell-annotation file.",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290578|https://pmc.ncbi.nlm.nih.gov/articles/PMC13267352/",
        source_accession="GSE290578;PRJNA1228650;phs001886.v1.p1", source_doi="10.1186/s13040-026-00556-y",
    ),
    pe_row(
        study_id="GSE298602", title="VGLL3-centered network connects placental, vascular, and immune defects in preeclampsia",
        year=2025, publication_status="PREPRINT_AS_OF_CUTOFF", pmid="40502186", doi="10.1101/2025.05.30.657097",
        pe_subtype="4/8 severe features; remaining subtype not specified", pe_donor_n=8, control_donor_n=3,
        tissue="placenta with chorioamniotic membranes", compartment="whole-thickness placenta and membranes",
        ga_case="mean 267 days, range 251-284", ga_control="mean 277.3 days, range 273-285", ga_matching="PARTIAL; small control group and case range includes preterm",
        platform="10x Chromium 3' v3; NovaSeq 6000", repository_accession="GSE298602; PRJNA1270014",
        processed_data_availability="PUBLIC_PROCESSED: raw and filtered matrices per GSM; no stand-alone annotation file found",
        raw_data_availability="PUBLIC_RAW: SRA", donor_identity_availability="YES: 11 GSM/sample IDs", patient_level_pseudobulk_feasible="YES",
        possible_overlap="NONE_IDENTIFIED", recommended_role="SCRNA_REPLICATION",
        evidence_summary="Paper calls 8 donors PreE, whereas GEO treatment labels split them as 4 PreE_SF and 4 gHTN; disease-label reconciliation is mandatory before reuse. Three controls limit discovery value.",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298602|https://pmc.ncbi.nlm.nih.gov/articles/PMC12157399/",
        source_accession="GSE298602;PRJNA1270014", source_pmid="40502186", source_doi="10.1101/2025.05.30.657097",
    ),
    pe_row(
        study_id="Yang2023_GitHub", title="Single-cell RNA-seq reveals developmental deficiencies in both placentation and decidualization in women with late-onset preeclampsia",
        year=2023, publication_status="PEER_REVIEWED", pmid="37283740", doi="10.3389/fimmu.2023.1142273",
        pe_subtype="LOPE", pe_donor_n="placenta 3; decidua 3", control_donor_n="placenta 3; decidua 4",
        tissue="placenta and decidua", compartment="near cord insertion; placenta and decidua analyzed separately",
        ga_case="placenta 37+2 to 37+6; decidua 35+1 to 37+6", ga_control="placenta 39+2 to 41+1; decidua 38+1 to 41+1",
        ga_matching="NO: LOPE earlier than controls", platform="custom Drop-seq; NovaSeq 6000",
        repository_accession="GitHub JustMoveOnnn/preeclampsia", processed_data_availability="PUBLIC_PROCESSED DGE matrices",
        raw_data_availability="NOT_FOUND", donor_identity_availability="PARTIAL: tissue sample labels and clinical table; no global donor key",
        patient_level_pseudobulk_feasible="YES WITHIN TISSUE; cross-tissue pairing only where clinical aliases match",
        possible_overlap="Clinical table supports deci3-pla2, deci5-pla3 and deci7-pla5 as paired pregnancies; other pairings are not supported",
        recommended_role="SCRNA_REPLICATION", evidence_summary="Direct spreadsheet inspection confirmed 6 placenta and 7 decidua biological samples. Repository file count must not be substituted for donor n. Suitable only for LOPE replication with tissue stratification and GA caveat.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC10239844/|https://github.com/JustMoveOnnn/preeclampsia/tree/main/single_cell_matrix/data",
        source_accession="GitHub:JustMoveOnnn/preeclampsia", source_pmid="37283740", source_doi="10.3389/fimmu.2023.1142273",
    ),
    pe_row(
        study_id="Tsang2017_EGAS00001002449", title="Integrative single-cell and cell-free plasma RNA transcriptomics elucidates placental cellular dynamics",
        year=2017, publication_status="PEER_REVIEWED", doi="10.1073/pnas.1710470114",
        pe_subtype="early/severe PE", pe_donor_n=4, control_donor_n=4, tissue="placenta", compartment="placental biopsies",
        ga_case="24 to 33+6 weeks", ga_control="normal term", ga_matching="NO: early PE vs term control",
        platform="10x Genomics", repository_accession="EGAS00001002449", processed_data_availability="PUBLIC_PROCESSED STATUS UNRESOLVED",
        raw_data_availability="CONTROLLED_RAW: EGA", donor_identity_availability="YES UNDER CONTROLLED ACCESS", patient_level_pseudobulk_feasible="CONDITIONAL",
        recommended_role="PENDING", evidence_summary="Foundational small cohort; controlled access and complete GA confounding prevent primary public discovery use.",
        source_url="https://doi.org/10.1073/pnas.1710470114|https://ega-archive.org/studies/EGAS00001002449",
        source_accession="EGAS00001002449", source_doi="10.1073/pnas.1710470114",
    ),
    pe_row(
        study_id="HRA003297", title="A Galectin-9-Driven CD11chigh Decidual Macrophage Subset Suppresses Uterine Vascular Remodeling in Preeclampsia",
        year=2024, publication_status="PEER_REVIEWED", pmid="38314577", doi="10.1161/CIRCULATIONAHA.123.064391",
        pe_donor_n=5, control_donor_n=3, tissue="decidua", compartment="maternal decidua; scRNA-seq plus spatial transcriptomics",
        platform="scRNA-seq and spatial transcriptomics", repository_accession="HRA003297; PRJCA012601",
        processed_data_availability="NOT_FOUND PUBLICLY", raw_data_availability="CONTROLLED_RAW",
        donor_identity_availability="BioSample aliases public; essential files controlled", patient_level_pseudobulk_feasible="CONDITIONAL ON ACCESS",
        recommended_role="PENDING", evidence_summary="NGDC lists eight BioSamples (5 PE, 3 normal) and Controlled Access. It cannot rescue public discovery without approved data access and a donor-clinical map.",
        source_url="https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003297|https://pubmed.ncbi.nlm.nih.gov/38314577/",
        source_accession="HRA003297;PRJCA012601", source_pmid="38314577", source_doi="10.1161/CIRCULATIONAHA.123.064391",
    ),
    pe_row(
        study_id="HRA004699", title="Single-cell profiling reveals immune disturbances landscape and HLA-F-mediated immune tolerance at the maternal-fetal interface in preeclampsia",
        year=2023, publication_status="PEER_REVIEWED", doi="10.3389/fimmu.2023.1234577",
        pe_donor_n=2, control_donor_n=2, tissue="decidua", compartment="maternal-fetal interface",
        ga_case="third trimester", ga_control="third trimester", ga_matching="INSUFFICIENT DETAIL",
        platform="BD Rhapsody; Illumina HiSeq 4000", repository_accession="HRA004699; PRJCA017294",
        processed_data_availability="NOT_FOUND PUBLICLY", raw_data_availability="CONTROLLED_RAW",
        donor_identity_availability="BioSample aliases public; files controlled", patient_level_pseudobulk_feasible="CONDITIONAL ON ACCESS",
        recommended_role="PENDING", evidence_summary="Only 2+2 donors and controlled access; appropriate at most as targeted decidual replication after access.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC10579943/|https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004699",
        source_accession="HRA004699;PRJCA017294", source_doi="10.3389/fimmu.2023.1234577",
    ),
    pe_row(
        study_id="GSE173193", title="Human placenta single-cell transcriptome cohort including PE, GDM and advanced maternal age",
        year=2021, publication_status="PEER_REVIEWED COMPONENT REPORTS", pe_subtype="severe PE", pe_donor_n=2, control_donor_n=2,
        other_donor_n="GDM 2; advanced maternal age 2", tissue="placenta", compartment="placental tissue",
        ga_case="32-40 weeks across cohort; PE-specific detail limited", ga_control="32-40 weeks across cohort", ga_matching="SMALL N; not adequate for discovery",
        platform="10x Genomics", repository_accession="GSE173193", processed_data_availability="PUBLIC_PROCESSED", raw_data_availability="PUBLIC_RAW",
        donor_identity_availability="YES", patient_level_pseudobulk_feasible="YES BUT n=2+2", recommended_role="SCRNA_REPLICATION",
        evidence_summary="Exactly 2 control, 2 PE, 2 GDM and 2 advanced-age donors. Cell counts do not enlarge biological n.",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173193", source_accession="GSE173193",
    ),
    pe_row(
        study_id="GSE282038_cross_series", title="Single-cell RNA sequencing of placenta in early-onset preeclampsia",
        year=2025, publication_status="PEER_REVIEWED", pe_subtype="EOPE with FGR", pe_donor_n="GSE282038 n=1; later GSE298119 adds 2",
        control_donor_n="0 internal; GSE267340 supplies 1; GSE298119 adds 1", tissue="placenta", compartment="placental tissue",
        ga_case="index 31+4 weeks", ga_control="external index 38+4 weeks", ga_matching="NO; disease, GA and series are confounded",
        platform="10x Genomics 3' v3.1", repository_accession="GSE282038; GSE267340; GSE298119",
        processed_data_availability="PUBLIC_PROCESSED", raw_data_availability="PUBLIC_RAW", donor_identity_availability="YES BY SERIES",
        patient_level_pseudobulk_feasible="YES TECHNICALLY; INVALID AS PRIMARY CONTRAST", possible_overlap="Controls/cases are split across three series",
        recommended_role="SCRNA_REPLICATION", evidence_summary="No internal control and EOPE is inseparable from FGR and cross-series/GA effects.",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282038|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267340|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298119",
        source_accession="GSE282038;GSE267340;GSE298119",
    ),
    pe_row(
        study_id="GSE329173", title="Single-cell transcriptomes of severe preeclampsia placentas",
        year=2026, publication_status="GEO_RECORD", pe_subtype="severe PE", pe_donor_n=3, control_donor_n=0,
        tissue="placenta", compartment="placental tissue", ga_matching="NO INTERNAL CONTROL", repository_accession="GSE329173",
        processed_data_availability="PUBLIC_PROCESSED", raw_data_availability="PUBLIC_RAW", donor_identity_availability="YES",
        patient_level_pseudobulk_feasible="NO CASE-CONTROL CONTRAST", recommended_role="SCRNA_REPLICATION",
        evidence_summary="Exactly three severe-PE donors and no internal control: external PE cell-state validation only.",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329173", source_accession="GSE329173",
    ),
    pe_row(
        study_id="Tsuda2024_HUM0443", title="Decidual CD4-positive T-cell single-cell profiling in early-onset preeclampsia",
        year=2024, publication_status="PEER_REVIEWED", doi="10.3389/fimmu.2024.1401738", pe_subtype="EOPE", pe_donor_n=3,
        control_donor_n="healthy early n=4; healthy late n=3", tissue="decidua", compartment="sorted CD4+ T cells only",
        ga_case="diagnosed <34 weeks; delivered 34-37 weeks", ga_control="early gestation and healthy term groups", ga_matching="NO single ideal matched control group",
        platform="single-cell RNA-seq", repository_accession="HUM0443.v1; DRA017833; E-GEAD-674",
        processed_data_availability="PUBLIC/REGISTERED PROCESSED COUNTS", raw_data_availability="PUBLIC/REGISTERED RAW",
        donor_identity_availability="YES", patient_level_pseudobulk_feasible="YES FOR SORTED CD4 CELLS", recommended_role="SCRNA_REPLICATION",
        evidence_summary="Useful only for a restricted decidual CD4 state; not a whole-tissue discovery cohort.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC11106458/|https://humandbs.dbcls.jp/en/hum0443-v1",
        source_accession="HUM0443.v1;DRA017833;E-GEAD-674", source_doi="10.3389/fimmu.2024.1401738",
    ),
    pe_row(
        study_id="Jiao2023", title="Dissecting human placental cell heterogeneity in preeclampsia and gestational diabetes using single-cell sequencing",
        year=2023, publication_status="PEER_REVIEWED", doi="10.1016/j.molimm.2023.07.005", pe_donor_n=5, control_donor_n=5, other_donor_n="GDM n=5",
        tissue="placenta", compartment="placental tissue", platform="scRNA-seq; 96,048 reported cells", repository_accession="NOT_FOUND",
        processed_data_availability="NOT_FOUND", raw_data_availability="NOT_FOUND", donor_identity_availability="NO PUBLIC DATA",
        patient_level_pseudobulk_feasible="NO", recommended_role="EXCLUDE", exclusion_reason="EXCLUDE_DATA_UNAVAILABLE",
        evidence_summary="Focused GEO/SRA/CNSA/NGDC/Figshare/GitHub searches found no reusable repository; publisher endpoint exposed bibliographic metadata only. Publication-only results are not reusable data.",
        source_url="https://doi.org/10.1016/j.molimm.2023.07.005", source_doi="10.1016/j.molimm.2023.07.005",
    ),
    pe_row(
        study_id="Zhang2021", title="Dissecting human trophoblast cell transcriptional heterogeneity in preeclampsia using single-cell RNA sequencing",
        year=2021, publication_status="PEER_REVIEWED", pmid="34212522", doi="10.1002/mgg3.1730", pe_donor_n=3, control_donor_n=3,
        tissue="placental parenchyma", compartment="2 cm deep and 5 cm from cord insertion",
        ga_case="34+5, 35+3, 35+1", ga_control="38, 38+2, 38+5", ga_matching="NO",
        platform="Illumina HiSeq X", repository_accession="AVAILABLE_ON_REQUEST ONLY",
        processed_data_availability="NOT_FOUND: Data S1-S8 are cell counts, DEG lists and enrichment tables, not matrices",
        raw_data_availability="NOT PUBLIC; available on request", donor_identity_availability="NO IN EXPRESSION OBJECT: aggregate GT01/XG01 labels",
        patient_level_pseudobulk_feasible="NO", recommended_role="EXCLUDE", exclusion_reason="EXCLUDE_DATA_UNAVAILABLE",
        evidence_summary="Supplementary files were directly inspected; none is a reusable count/expression matrix. 'Available on request' is not public reuse.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC8404237/", source_pmid="34212522", source_doi="10.1002/mgg3.1730",
    ),
    pe_row(
        study_id="Nonn2025_HypertensionResearch", title="Single-nucleus transcriptomics in preeclampsia",
        year=2025, publication_status="PEER_REVIEWED", pmid="41419624", doi="10.1038/s41440-025-02504-5", pe_subtype="EOPE and LOPE",
        pe_donor_n=4, control_donor_n=4, tissue="placenta", compartment="villi and decidua",
        repository_accession="AVAILABLE_ON_REQUEST ONLY", processed_data_availability="NOT PUBLIC", raw_data_availability="NOT PUBLIC",
        donor_identity_availability="NO PUBLIC DATA", patient_level_pseudobulk_feasible="NO", recommended_role="EXCLUDE", exclusion_reason="EXCLUDE_DATA_UNAVAILABLE",
        evidence_summary="Eight samples (2 EOPE, 2 LOPE, 2 preterm control, 2 term control) are scientifically relevant but not publicly downloadable.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12960252/", source_pmid="41419624", source_doi="10.1038/s41440-025-02504-5",
    ),
    pe_row(
        study_id="Mu2026", title="Revealing an enhanced cytotoxic immune microenvironment at the human maternal-fetal interface in preeclampsia",
        year=2026, publication_status="PEER_REVIEWED", pmid="41991681", doi="10.1007/s10565-026-10170-7", pe_subtype="primarily LOPE",
        pe_donor_n=3, control_donor_n=3, tissue="decidua", compartment="sorted CD45+ immune cells",
        repository_accession="NOT_FOUND", processed_data_availability="NOT_FOUND", raw_data_availability="NOT_FOUND",
        donor_identity_availability="NO PUBLIC DATA", patient_level_pseudobulk_feasible="NO", recommended_role="EXCLUDE", exclusion_reason="EXCLUDE_DATA_UNAVAILABLE",
        evidence_summary="The paper reports 3+3 sorted decidual immune samples but no reusable public accession was located by the cutoff.",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC13315417/", source_pmid="41991681", source_doi="10.1007/s10565-026-10170-7",
    ),
    pe_row(
        study_id="Wei2025", title="Single-cell RNA sequencing reveals systemic and placental immune landscape in preeclampsia",
        year=2025, publication_status="PEER_REVIEWED", doi="10.1016/j.placenta.2025.08.325", pe_subtype="EOPE and LOPE",
        pe_donor_n="EOPE 2; LOPE 2", control_donor_n=4, tissue="placenta and peripheral blood", compartment="multiple placental compartments and PBMC",
        repository_accession="NOT_FOUND", processed_data_availability="NOT_FOUND", raw_data_availability="NOT_FOUND",
        donor_identity_availability="NO PUBLIC DATA", patient_level_pseudobulk_feasible="NO", recommended_role="EXCLUDE", exclusion_reason="EXCLUDE_DATA_UNAVAILABLE",
        evidence_summary="No reusable repository was located; publication-only findings are excluded from dataset counts.",
        source_url="https://doi.org/10.1016/j.placenta.2025.08.325", source_doi="10.1016/j.placenta.2025.08.325",
    ),
    pe_row(
        study_id="Xiao2025", title="Maternal-Fetal Interface Cell Dysfunction in Patients With Preeclampsia Revealed via Single-Cell RNA Sequencing",
        year=2025, publication_status="PEER_REVIEWED", doi="10.1111/aji.70101", pe_donor_n=3, control_donor_n=3,
        tissue="placenta and decidua", compartment="maternal-fetal interface tissue after cesarean delivery",
        platform="single-cell RNA-seq; 32,279 reported cells", repository_accession="NOT_FOUND",
        processed_data_availability="NOT_FOUND", raw_data_availability="NOT_FOUND", donor_identity_availability="NO PUBLIC DATA",
        patient_level_pseudobulk_feasible="NO", recommended_role="EXCLUDE", exclusion_reason="EXCLUDE_DATA_UNAVAILABLE",
        evidence_summary="The paper reports 3 PE+3 normal pregnancies, but no reusable public repository or expression matrix was located by the cutoff.",
        source_url="https://doi.org/10.1111/aji.70101", source_doi="10.1111/aji.70101",
    ),
    pe_row(
        study_id="GSE265862", title="Mapping Decidualization Resistance in Former Severe Preeclampsia Patients at Multi-Omic Levels",
        year=2024, publication_status="PEER_REVIEWED", pe_subtype="history of severe PE", pe_donor_n=11, control_donor_n=12,
        tissue="nonpregnant endometrium", compartment="late-secretory endometrial biopsy after prior pregnancy",
        repository_accession="GSE265862", processed_data_availability="PUBLIC_PROCESSED", raw_data_availability="PUBLIC_RAW",
        donor_identity_availability="YES", patient_level_pseudobulk_feasible="TECHNICALLY YES BUT OUT OF SCOPE", recommended_role="EXCLUDE", exclusion_reason="OUT_OF_SCOPE_TISSUE_AND_TIMEPOINT",
        evidence_summary="Not a current PE placenta or maternal-fetal-interface case-control cohort; retained only as a recall-check exclusion.",
        source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265862", source_accession="GSE265862",
    ),
]


ACCESS_ROWS = [
    {"study_id":"ADMati_2023_FIGSHARE","component":"all-cell expression plus metadata/annotations","access_class":"PUBLIC_PROCESSED","essential_for_reanalysis":"YES","authorized_access_required":"NO","downloadable_now":"YES","patient_ids_retained":"YES","details":"256,501,211-byte ZIP; 5.87 GB uncompressed; 86,752 cells; 26 donorID values","source_url":"https://doi.org/10.6084/m9.figshare.23264102.v1","source_accession":"Figshare 23264102.v1","audit_date":AUDIT_DATE},
    {"study_id":"ADMati_2023_FIGSHARE","component":"raw FASTQ","access_class":"NOT_FOUND","essential_for_reanalysis":"NO for processed-matrix pseudobulk","authorized_access_required":"UNRESOLVED","downloadable_now":"NO","patient_ids_retained":"NOT_APPLICABLE","details":"No public raw accession located","source_url":"https://doi.org/10.1016/j.medj.2023.07.005","source_accession":"NOT_FOUND","audit_date":AUDIT_DATE},
    {"study_id":"PMID41472684_phs001886v6","component":"exact GEO processed/count matrix","access_class":"NOT_FOUND","essential_for_reanalysis":"YES unless controlled raw is obtained","authorized_access_required":"NO PUBLIC OBJECT EXISTS","downloadable_now":"NO","patient_ids_retained":"NO","details":"Exact-title GEO query returned zero; manuscript provides no GSE accession","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12747280/","source_accession":"GEO=NOT_FOUND","audit_date":AUDIT_DATE},
    {"study_id":"PMID41472684_phs001886v6","component":"raw sequence for the new 78-person cohort","access_class":"CONTROLLED_RAW","essential_for_reanalysis":"YES","authorized_access_required":"YES","downloadable_now":"NOT CONFIRMED: v6 expected sequence count is not represented in public SSTR records","patient_ids_retained":"NO PUBLIC CASE MAP","details":"dbGaP authorization would be required; actual v6 sequence loading/release remains unresolved","source_url":"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001886.v6.p1","source_accession":"phs001886.v6.p1","audit_date":AUDIT_DATE},
    {"study_id":"PMID41472684_phs001886v6","component":"individual clinical phenotypes and maternal/fetal mappings","access_class":"CONTROLLED_CLINICAL","essential_for_reanalysis":"YES","authorized_access_required":"YES","downloadable_now":"NO PUBLICLY","patient_ids_retained":"PUBLICLY NO","details":"Only aggregate Supplementary Table 2 is public","source_url":"https://doi.org/10.21203/rs.3.rs-8254581/v1","source_accession":"phs001886.v6.p1","audit_date":AUDIT_DATE},
    {"study_id":"GSE290578","component":"Cell Ranger matrices","access_class":"PUBLIC_PROCESSED","essential_for_reanalysis":"YES","authorized_access_required":"NO","downloadable_now":"YES","patient_ids_retained":"YES","details":"Per-GSM barcodes/features/matrix; no stand-alone annotation file","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290578","source_accession":"GSE290578","audit_date":AUDIT_DATE},
    {"study_id":"GSE290578","component":"FASTQ","access_class":"PUBLIC_RAW","essential_for_reanalysis":"NO","authorized_access_required":"NO","downloadable_now":"YES","patient_ids_retained":"YES","details":"16 SRA-linked tissue libraries from 8 pregnancies","source_url":"https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1228650","source_accession":"PRJNA1228650","audit_date":AUDIT_DATE},
    {"study_id":"GSE298602","component":"raw/filtered matrices and FASTQ","access_class":"PUBLIC_PROCESSED+PUBLIC_RAW","essential_for_reanalysis":"YES","authorized_access_required":"NO","downloadable_now":"YES","patient_ids_retained":"YES","details":"11 GSM; 4 PreE_SF, 4 gHTN labels and 3 controls","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298602","source_accession":"GSE298602;PRJNA1270014","audit_date":AUDIT_DATE},
    {"study_id":"Yang2023_GitHub","component":"DGE matrices and clinical table","access_class":"PUBLIC_PROCESSED","essential_for_reanalysis":"YES","authorized_access_required":"NO","downloadable_now":"YES","patient_ids_retained":"PARTIAL","details":"Tissue-specific labels; clinical matching supports only documented cross-tissue pairs","source_url":"https://github.com/JustMoveOnnn/preeclampsia/tree/main/single_cell_matrix/data","source_accession":"GitHub repository","audit_date":AUDIT_DATE},
    {"study_id":"HRA003297","component":"raw scRNA/spatial files","access_class":"CONTROLLED_RAW","essential_for_reanalysis":"YES","authorized_access_required":"YES","downloadable_now":"NO WITHOUT APPROVAL","patient_ids_retained":"CONDITIONAL","details":"5 PE and 3 normal BioSamples","source_url":"https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003297","source_accession":"HRA003297","audit_date":AUDIT_DATE},
    {"study_id":"HRA004699","component":"raw scRNA files","access_class":"CONTROLLED_RAW","essential_for_reanalysis":"YES","authorized_access_required":"YES","downloadable_now":"NO WITHOUT APPROVAL","patient_ids_retained":"CONDITIONAL","details":"2 PE and 2 normal donors","source_url":"https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004699","source_accession":"HRA004699","audit_date":AUDIT_DATE},
    {"study_id":"Jiao2023","component":"expression matrix/raw data","access_class":"NOT_FOUND","essential_for_reanalysis":"YES","authorized_access_required":"UNKNOWN","downloadable_now":"NO","patient_ids_retained":"NO","details":"No reusable repository located","source_url":"https://doi.org/10.1016/j.molimm.2023.07.005","source_accession":"NOT_FOUND","audit_date":AUDIT_DATE},
    {"study_id":"Zhang2021","component":"expression matrix/raw data","access_class":"NOT_FOUND","essential_for_reanalysis":"YES","authorized_access_required":"NO; request to authors would be required","downloadable_now":"NO","patient_ids_retained":"NO","details":"Supplements contain summary tables only","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC8404237/","source_accession":"AVAILABLE_ON_REQUEST","audit_date":AUDIT_DATE},
]


HUCMSC_ROWS = [
    {"dataset":"GSE182158","independent_donor_n":"2 UC donors (U01,U02); 11 MSC donors total","tissue_definition":"umbilical-cord MSC; Wharton's jelly not explicitly resolved","passage":"P1 or P2 protocol; donor-specific passage unresolved","culture_condition":"cultured, unstimulated","inflammatory_stimulation":"NONE","donor_sex_age":"U01 female 28; U02 female 37","platform":"10x Chromium 3' v2; NovaSeq 6000","raw_availability":"PUBLIC_RAW","processed_availability":"PUBLIC_PROCESSED","donor_identity":"YES","sender_use":"BASELINE CORE, but never alone","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182158|https://pmc.ncbi.nlm.nih.gov/articles/PMC8715893/","source_accession":"GSE182158","audit_date":AUDIT_DATE},
    {"dataset":"GSE199071","independent_donor_n":"4 HUCMSC donors","tissue_definition":"human umbilical-cord MSC; exact WJ/whole-cord per donor unresolved","passage":"P3 and P6 labels","culture_condition":"cultured, unstimulated HUCMSC; HUVEC excluded","inflammatory_stimulation":"NONE","donor_sex_age":"HUCMSC1 female; HUCMSC2-4 male; age unresolved","platform":"10x scRNA-seq; NextSeq 550","raw_availability":"PUBLIC_RAW","processed_availability":"PUBLIC_PROCESSED","donor_identity":"YES","sender_use":"INDEPENDENT BASELINE REPLICATION","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199071","source_accession":"GSE199071","audit_date":AUDIT_DATE},
    {"dataset":"GSE117837_UNTREATED","independent_donor_n":"2 UC donors","tissue_definition":"umbilical-cord-derived MSC; exact WJ/whole-cord unresolved","passage":"D1 P5; D2 P0/P2/P5","culture_condition":"cultured untreated strata only","inflammatory_stimulation":"NONE in baseline arm","donor_sex_age":"UNRESOLVED","platform":"Fluidigm C1; HiSeq 2500","raw_availability":"PUBLIC_RAW","processed_availability":"PUBLIC_PROCESSED","donor_identity":"YES","sender_use":"PASSAGE-SENSITIVITY AND BASELINE REPLICATION","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837|https://pmc.ncbi.nlm.nih.gov/articles/PMC6506509/","source_accession":"GSE117837","audit_date":AUDIT_DATE},
    {"dataset":"GSE117837_LICENSED","independent_donor_n":"2 UC donors; valid paired strata D1 P5 and D2 P2/P5","tissue_definition":"umbilical-cord-derived MSC","passage":"D1 P5; D2 P2/P5; D2 P0 has no stimulated match","culture_condition":"same-donor/same-passage paired strata","inflammatory_stimulation":"IFN-gamma 10 ng/mL + TNF-alpha 10 ng/mL for 12 h","donor_sex_age":"UNRESOLVED","platform":"Fluidigm C1; HiSeq 2500","raw_availability":"PUBLIC_RAW","processed_availability":"PUBLIC_PROCESSED","donor_identity":"YES","sender_use":"LICENSING DESIGN ONLY; do not pool cells","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837|https://pmc.ncbi.nlm.nih.gov/articles/PMC6506509/","source_accession":"GSE117837","audit_date":AUDIT_DATE},
    {"dataset":"CNP0000562","independent_donor_n":"3 umbilical cords (UC1-UC3)","tissue_definition":"Wharton's jelly MSC from full-term natural births","passage":"expanded twice; exact conventional passage label unresolved","culture_condition":"cultured, unstimulated","inflammatory_stimulation":"NONE","donor_sex_age":"newborn sex 2 female/1 male; maternal age unresolved","platform":"10x scRNA-seq; 6,878 filtered cells","raw_availability":"CNGB project registered; direct file accessibility unresolved","processed_availability":"UNRESOLVED","donor_identity":"YES IN PAPER (UC1-UC3)","sender_use":"PENDING BASELINE REPLICATION UNTIL DOWNLOAD VERIFIED","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC7132901/|https://db.cngb.org/search/project/CNP0000562/","source_accession":"CNP0000562","audit_date":AUDIT_DATE},
    {"dataset":"HRA005090","independent_donor_n":"1 scRNA-seq donor (seven condition samples are not seven donors)","tissue_definition":"Wharton's jelly MSC","passage":"P4","culture_condition":"unprimed plus six cytokine-primed arms","inflammatory_stimulation":"24 h: IFNg 25 ng/mL; TNFa 10; IL4 50; IL6 100; IL15 100; IL17 100","donor_sex_age":"UNRESOLVED","platform":"10x scRNA-seq","raw_availability":"OPEN ACCESS GSA-Human","processed_availability":"UNRESOLVED","donor_identity":"YES, single donor","sender_use":"LICENSING SUPPLEMENT ONLY; cannot establish donor robustness","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC10964690/|https://ngdc.cncb.ac.cn/gsa-human/browse/HRA005090","source_accession":"HRA005090;PRJCA018383","audit_date":AUDIT_DATE},
    {"dataset":"PRJNA643879","independent_donor_n":"1 cord","tissue_definition":"uncultured Wharton's jelly/umbilical cord mixed cells","passage":"P0/uncultured","culture_condition":"fresh tissue cell suspension","inflammatory_stimulation":"NONE","donor_sex_age":"UNRESOLVED","platform":"10x scRNA-seq; 5,330 cells","raw_availability":"PUBLIC_RAW","processed_availability":"PUBLICATION-LEVEL; reusable matrix unresolved","donor_identity":"YES, single donor","sender_use":"TISSUE-ORIGIN REFERENCE ONLY; not cultured therapeutic sender","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC7791785/|https://www.ncbi.nlm.nih.gov/bioproject/PRJNA643879","source_accession":"PRJNA643879","audit_date":AUDIT_DATE},
    {"dataset":"WJMSC_bulk_2025","independent_donor_n":"7 donors","tissue_definition":"Wharton's jelly MSC from full-term spontaneous births","passage":"P3-P5","culture_condition":"naive vs licensed bulk RNA-seq","inflammatory_stimulation":"TNFa 10 ng/mL + IFNg 10 ng/mL for 24 h","donor_sex_age":"UNRESOLVED","platform":"bulk RNA-seq","raw_availability":"NOT_FOUND","processed_availability":"NOT_FOUND","donor_identity":"reported in paper only","sender_use":"METHOD NEAR-NEIGHBOR; EXCLUDE_DATA_UNAVAILABLE","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12010610/","source_accession":"NOT_FOUND","audit_date":AUDIT_DATE},
]


RISK_ROWS = [
    {"risk_id":"P0B-001","scope":"PMID41472684/phs001886.v6","severity":"CRITICAL","phase1_blocking":"YES for using this cohort","risk":"Exact GEO accession and public processed matrix not found","mitigation":"Keep PENDING; obtain exact GSE or dbGaP authorization plus manifest","evidence":"Exact-title GDS search count=0; preprint names only phs001886.v6","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12747280/","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-002","scope":"PMID41472684/phs001886.v6","severity":"CRITICAL","phase1_blocking":"YES for using this cohort","risk":"Public SSTR cannot map 78 patients to case/subtype/clinical fields","mitigation":"Require subject-sample-clinical crosswalk before pseudobulk","evidence":"v6 mixes maternal/fetal genotype and other records; 374 subjects/395 samples","source_url":"https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001886.v6.p1","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-003","scope":"phs001886.v1-v6","severity":"HIGH","phase1_blocking":"YES for independence claims","risk":"All versions are cumulative, not independent cohorts","mitigation":"Treat versions as one evolving study family","evidence":"Exact public dbGaP IDs from each earlier release are nested in v6","source_url":"https://www.ncbi.nlm.nih.gov/gap/sstr/","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-004","scope":"GSE290578","severity":"HIGH","phase1_blocking":"NO if replication only","risk":"Complete disease/gestational-age separation: EOPE preterm vs controls term","mitigation":"Replication/cell-state concordance only; never causal disease discovery","evidence":"GEO/paper GA 28+5-34+1 vs 38+2-39+6","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC13267352/","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-005","scope":"GSE290578","severity":"HIGH","phase1_blocking":"YES for naive n=16 analysis","risk":"Two placental layers are paired within each of eight pregnancies","mitigation":"Pregnancy is the replicate; layer is paired/repeated tissue","evidence":"16 GSM encode MF/Pla pairs for Norm1/4/6/7 and PePT2/3/5/8","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290578","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-006","scope":"GSE290578","severity":"HIGH","phase1_blocking":"NO","risk":"Near-neighbor already performs PE scRNA, ML, ligand/receptor, cell communication and SPP1-CD44 prioritization","mitigation":"Novelty must be therapeutic sender-to-PE target validation, not generic communication prioritization","evidence":"Peer-reviewed BioData Mining methods/results","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC13267352/","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-007","scope":"ADMati_2023_FIGSHARE","severity":"MODERATE","phase1_blocking":"NO","risk":"31 libraries represent 26 donors and include EOPE/LOPE/IUGR heterogeneity","mitigation":"Collapse by donorID; stratify subtype/IUGR; model GA","evidence":"Direct processed-matrix header inspection","source_url":"https://doi.org/10.6084/m9.figshare.23264102.v1","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-008","scope":"GSE298602","severity":"HIGH","phase1_blocking":"YES until labels reconciled","risk":"Paper calls 8 donors PreE but GEO labels 4 PreE_SF and 4 gHTN","mitigation":"Resolve clinical diagnosis mapping before PE replication","evidence":"GEO SOFT vs preprint Table 1","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298602","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-009","scope":"Yang2023","severity":"HIGH","phase1_blocking":"NO if tissue-specific replication","risk":"Placenta/decidua contain paired pregnancies and strong GA mismatch","mitigation":"Never add tissues as independent donors; analyze tissues separately","evidence":"Direct Table S1 clinical-profile matching","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC10239844/","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-010","scope":"GSE182158","severity":"HIGH","phase1_blocking":"YES for two-donor-only sender","risk":"UC arm contains only U01 and U02","mitigation":"Require cross-dataset donor-level replication","evidence":"GEO/paper Supplementary Table 1","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182158","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-011","scope":"GSE117837","severity":"HIGH","phase1_blocking":"YES for pooled cell comparison","risk":"Donor and passage strata are unbalanced; D2 P0 has no stimulated match","mitigation":"Only same-donor/same-passage contrasts D1 P5, D2 P2 and D2 P5","evidence":"GEO per-cell labels","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-012","scope":"GSE117837","severity":"MODERATE","phase1_blocking":"YES for exact cell-count claims","risk":"GEO labels 203 naive/158 stimulated; paper states 202/159","mitigation":"Retain submitted labels and resolve one-cell discrepancy before analysis","evidence":"GEO vs publication","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC6506509/","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-013","scope":"GSE272342","severity":"CRITICAL","phase1_blocking":"YES for general replication claim","risk":"All controls are twin placentas; cases mix singleton and twin pregnancies","mitigation":"Use only design-specific subtype/sensitivity analyses at pregnancy level","evidence":"GEO overall design","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272342","audit_date":AUDIT_DATE},
    {"risk_id":"P0B-014","scope":"Jiao2023/Zhang2021/other publication-only cohorts","severity":"HIGH","phase1_blocking":"YES for inclusion","risk":"Publication results are not reusable expression data","mitigation":"EXCLUDE_DATA_UNAVAILABLE; do not count toward replication","evidence":"Repository search and supplement inspection","source_url":"https://doi.org/10.1016/j.molimm.2023.07.005|https://pmc.ncbi.nlm.nih.gov/articles/PMC8404237/","audit_date":AUDIT_DATE},
]


ROLE_ROWS = [
    {"dataset":"Admati_2023_FIGSHARE","phase0a_role":"NOT_AUDITED","phase0b_role":"PRIMARY_PE_SCRNA","change":"ADDED/RESCUES PRIMARY DISCOVERY","rationale":"Public 86,752-cell matrix retains 26 donor IDs and subtype/clinical fields","mandatory_restrictions":"Donor-level pseudobulk; subtype/IUGR/GA-aware; 31 libraries are not 31 donors","source_url":"https://doi.org/10.6084/m9.figshare.23264102.v1","audit_date":AUDIT_DATE},
    {"dataset":"PMID41472684_phs001886.v6","phase0a_role":"NOT_AUDITED","phase0b_role":"PENDING","change":"ADDED","rationale":"Reported 39+39 but essential public matrix/GEO and patient mapping not found","mandatory_restrictions":"Do not use until access and 78-subject crosswalk are independently verified","source_url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12747280/","audit_date":AUDIT_DATE},
    {"dataset":"GSE290578","phase0a_role":"NOT_AUDITED","phase0b_role":"SCRNA_REPLICATION","change":"ADDED","rationale":"Public paired layers from 4 PE+4 controls, but severe GA confounding","mandatory_restrictions":"Pregnancy-level replication; paired layer; not primary discovery","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290578","audit_date":AUDIT_DATE},
    {"dataset":"GSE298602","phase0a_role":"NOT_AUDITED","phase0b_role":"SCRNA_REPLICATION","change":"ADDED","rationale":"Public 8 reported PreE+3 controls; small controls and label inconsistency","mandatory_restrictions":"Resolve PreE_SF versus gHTN labels first","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298602","audit_date":AUDIT_DATE},
    {"dataset":"Yang2023_GitHub","phase0a_role":"NOT_AUDITED","phase0b_role":"SCRNA_REPLICATION","change":"ADDED","rationale":"Public LOPE processed matrices with tissue-level donor labels","mandatory_restrictions":"Tissue-stratified; paired pregnancies not independent; GA caveat","source_url":"https://github.com/JustMoveOnnn/preeclampsia/tree/main/single_cell_matrix/data","audit_date":AUDIT_DATE},
    {"dataset":"GSE173193","phase0a_role":"SCRNA_REPLICATION","phase0b_role":"SCRNA_REPLICATION","change":"UNCHANGED","rationale":"Only 2 PE+2 control donors","mandatory_restrictions":"Auxiliary validation only","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173193","audit_date":AUDIT_DATE},
    {"dataset":"GSE282038/GSE267340/GSE298119","phase0a_role":"SCRNA_REPLICATION","phase0b_role":"SCRNA_REPLICATION","change":"UNCHANGED","rationale":"Cross-series controls and EOPE-FGR/GA confounding","mandatory_restrictions":"No primary discovery contrast","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282038","audit_date":AUDIT_DATE},
    {"dataset":"GSE329173","phase0a_role":"SCRNA_REPLICATION","phase0b_role":"SCRNA_REPLICATION","change":"UNCHANGED","rationale":"Three severe PE samples; no internal control","mandatory_restrictions":"External PE cell-state validation only","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329173","audit_date":AUDIT_DATE},
    {"dataset":"GSE182158","phase0a_role":"PRIMARY_HUCMSC_ATLAS","phase0b_role":"PRIMARY_HUCMSC_ATLAS","change":"RESTRICTED","rationale":"Good atlas but only two UC donors","mandatory_restrictions":"Cannot define robust sender alone; require GSE199071 and GSE117837 donor-level replication","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182158","audit_date":AUDIT_DATE},
    {"dataset":"GSE199071","phase0a_role":"SUPPLEMENTARY","phase0b_role":"SUPPLEMENTARY","change":"PROMOTED WITHIN SENDER EVIDENCE","rationale":"Four independent HUCMSC donor labels provide baseline cross-dataset replication","mandatory_restrictions":"Account for P3/P6; exclude HUVEC","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199071","audit_date":AUDIT_DATE},
    {"dataset":"GSE117837","phase0a_role":"HUCMSC_LICENSING","phase0b_role":"HUCMSC_LICENSING","change":"RESTRICTED DESIGN","rationale":"Valid within-stratum contrasts exist only for D1 P5, D2 P2, D2 P5","mandatory_restrictions":"No pooled cell-level comparison; D2 P0 is baseline only","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117837","audit_date":AUDIT_DATE},
    {"dataset":"GSE272342","phase0a_role":"BULK_REPLICATION","phase0b_role":"SUBTYPE_VALIDATION","change":"DOWNGRADED","rationale":"All controls are twin placentas while cases mix singleton/twin; multiplicity is structurally confounded","mandatory_restrictions":"DESIGN_SPECIFIC_SENSITIVITY only; pregnancy-level modeling; preserve GSE203507 overlap flag","source_url":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272342","audit_date":AUDIT_DATE},
    {"dataset":"Nature2026_normal_atlas","phase0a_role":"NOT_AUDITED","phase0b_role":"SUPPLEMENTARY","change":"ADDED REFERENCE","rationale":"Normal MFI across GW5-39 with public annotated/raw-count h5ad","mandatory_restrictions":"NORMAL_GESTATIONAL_REFERENCE/CELL_ANNOTATION_REFERENCE only; never PE case-control","source_url":"https://doi.org/10.1038/s41586-026-10316-x|https://cell.ucsf.edu/snPlacenta/","audit_date":AUDIT_DATE},
]


def dbgap_sets(version: str) -> tuple[set[int], set[int], set[str], set[str]]:
    payload = json.loads((RAW / f"{version}_sstr_subjects.json").read_text(encoding="utf-8"))
    rows = payload["data"]
    subjects = {int(row["dbgap_subject_id"]) for row in rows if row.get("dbgap_subject_id") is not None}
    samples = {int(row["dbgap_sample_id"]) for row in rows if row.get("dbgap_sample_id") is not None}
    subject_aliases = {str(row["submitted_subject_id"]) for row in rows if row.get("submitted_subject_id")}
    sample_aliases = {str(row["submitted_sample_id"]) for row in rows if row.get("submitted_sample_id")}
    return subjects, samples, subject_aliases, sample_aliases


def build_overlap() -> list[dict[str, object]]:
    versions = [f"phs001886.v{i}.p1" for i in range(1, 7)]
    sets = {version: dbgap_sets(version) for version in versions}
    rows: list[dict[str, object]] = []
    for i, left in enumerate(versions):
        for right in versions[i + 1:]:
            ls, lm, lsa, lma = sets[left]
            rs, rm, rsa, rma = sets[right]
            rows.append({
                "earlier_version": left, "later_version": right,
                "earlier_subject_n": len(ls), "later_subject_n": len(rs), "shared_subject_n": len(ls & rs),
                "new_subject_n": len(rs - ls), "retired_subject_n": len(ls - rs),
                "earlier_sample_n": len(lm), "later_sample_n": len(rm), "shared_sample_n": len(lm & rm),
                "new_sample_n": len(rm - lm), "retired_sample_n": len(lm - rm),
                "shared_submitted_subject_alias_n": len(lsa & rsa), "shared_submitted_sample_alias_n": len(lma & rma),
                "independence_conclusion": "NOT_INDEPENDENT_CUMULATIVE" if ls <= rs and lm <= rm else "OVERLAPPING_NOT_STRICTLY_NESTED",
                "source_url": "https://www.ncbi.nlm.nih.gov/gap/sstr/", "source_accession": f"{left};{right}",
                "audit_date": AUDIT_DATE,
            })
    return rows


def revise_phase0a_gse272342() -> None:
    for relative in ["results/00_dataset_audit/dataset_registry.csv", "results/00_dataset_audit/proposed_dataset_roles.csv"]:
        path = ROOT / relative
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = list(reader)
        for row in rows:
            if row.get("geo_accession") == "GSE272342":
                row["proposed_role"] = "SUBTYPE_VALIDATION"
                if "role_restrictions" in row:
                    row["role_restrictions"] = "DESIGN_SPECIFIC_SENSITIVITY only: all controls are twin placentas while cases mix singleton/twin; model pregnancy as the independent unit; stratify multiplicity/chorionicity; preserve three GSE203507 overlap flags"
                if "mandatory_restrictions" in row:
                    row["mandatory_restrictions"] = "DESIGN_SPECIFIC_SENSITIVITY only; all controls are twins and cases mix singleton/twin; pregnancy-level modeling; preserve three GSE203507 overlaps"
                if "rationale" in row:
                    row["rationale"] = "All controls are twin placentas while cases contain singleton and twin pregnancies; a general PE replication claim is structurally confounded by multiplicity."
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def main() -> int:
    write_csv("pe_scrna_extended_registry.csv", PE_ROWS)
    write_csv("pe_scrna_data_access.csv", ACCESS_ROWS)
    write_csv("phs001886_version_overlap.csv", build_overlap())
    write_csv("hucmsc_sender_redundancy_registry.csv", HUCMSC_ROWS)
    write_csv("phase0b_risk_flags.csv", RISK_ROWS)
    write_csv("revised_dataset_roles.csv", ROLE_ROWS)
    revise_phase0a_gse272342()
    print(f"Phase 0B registries written to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
