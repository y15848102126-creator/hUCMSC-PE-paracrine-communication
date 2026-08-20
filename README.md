# hUCMSC–PE therapeutic communication — public repository candidate v1.1

## Release status

**HUMAN-REVIEWED RELEASE CANDIDATE — PREPARED FOR PRIVATE GITHUB INSPECTION; NOT YET PUBLICLY PUBLISHED.**

This sanitized candidate reconstructs the frozen computational workflow accompanying the manuscript *Independent mapping of placental receiver states and hUC-MSC sender programs constrains putative paracrine communication in preeclampsia*. It does not establish a therapeutic mechanism.

## Authoritative and superseded methods

> **AUTHORITATIVE RECEIVER METHOD:** Phase 2A.2 pregnancy-level continuous-expression analysis.
>
> **SUPERSEDED / HISTORICAL ONLY:** Phase 2A count-likelihood edgeR receiver analysis. Historical files carry explicit deprecation headers and are excluded from the default execution workflow.
>
> **WITHDRAWN / NON-AUTHORITATIVE:** the Phase 1A GSE30186 arbitrary cohort-minimum shift-log route. Formal outcome analysis used the frozen Phase 1A.1 normexp → quantile normalization → log2 matrix.

## Default reproducibility workflow

1. Review `config/`, `environment/`, the dataset-reference registry, and `THIRD_PARTY_RESOURCE_RETRIEVAL.md`.
2. Obtain public source data and database resources from their original providers; raw matrices and third-party dumps are not redistributed here.
3. Use the six frozen Phase 1A.1 bulk matrices and the authoritative Phase 2A.2 corrected receiver route.
4. Continue through Phase 2B receiver validation, Phase 3 sender analysis, Phase 4 integration, and external triangulation.
5. Run phase-specific validation scripts and reconcile outputs against `results/final_synthesis/`.

## Non-inferential receptor provenance correction

The Phase 4A receptor export previously duplicated pooled and subtype-specific eligibility rows and misclassified pooled `PE` summary rows into descriptive controls. Version 1.1 reconstructs only the five control-description fields from nine unique EARLY_CONTROL/LATE_CONTROL pregnancies. All PE receptor values, receptor-competence classes, downstream Phase 4A classifications, final candidate classifications, manuscript numbers, and conclusions are unchanged.

## Redistribution boundary

- No raw expression matrices, FASTQ files, controlled clinical data, source publication PDFs/supplements, complete third-party networks, or source proteomic identification tables are included.
- Complete MSigDB gene membership lists are omitted. The release retains gene-set identifiers, collection/version, sizes, membership hashes, source URLs, and retrieval/reconstruction instructions.
- Project-derived classifications, aggregate statistics, identifiers, hashes, and source/version records are retained for auditability.
- EV/exosome detection remains distinct from soluble conditioned-medium evidence.

## Licensing

The BSD 3-Clause License in `LICENSE` applies only to original project code. It does not relicense GEO, MSigDB, NicheNet, OmniPath, ProteomeXchange, publication-derived, or other third-party content. Derived tables remain subject to the boundaries described in `THIRD_PARTY_NOTICE.md` and `DATA_AND_DERIVED_OUTPUTS_NOTICE.md`; users must obtain source resources from their original providers and comply with provider terms.
