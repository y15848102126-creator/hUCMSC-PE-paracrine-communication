# Public GitHub release preparation report

Preparation date: 2026-08-20

Repository name: `hUCMSC-PE-paracrine-communication`
Recommended description: Reproducible code and derived outputs for cross-dataset inference of hUC-MSC sender–receiver communication in preeclampsia.

## Repository status

| Field | Status |
|---|---|
| Local repository path | `release/public_repository_candidate` within the private project workspace; the user-specific absolute path is intentionally omitted from the public report |
| GitHub repository URL | NOT CREATED — GitHub CLI/authentication unavailable during preparation |
| Visibility | LOCAL ONLY; no public release authorized or performed |
| Current branch | `main` |
| Initial release commit SHA | `d12bfd3cb356ef5aebcd3cb469465e5eb68e8280` |
| Local tag status | `v1.0.0` created locally after final repository verification; not pushed and not published as a GitHub Release |
| GitHub Release status | NOT CREATED; requires human approval after repository inspection |
| License | BSD 3-Clause for original project code only |
| Tracked files | 160, including this report and excluding `.git/` |
| Repository content size | 6,157,946 bytes, excluding `.git/` |
| Files excluded by `.gitignore` | 0 existing release-candidate files |
| Privacy/credential scan | PASS before Git initialization; no unresolved sensitive/private content detected |
| Manifest verification | PASS after regeneration; `FILE_MANIFEST.csv` excludes itself and covers every other tracked release file |
| Third-party redistribution status | No raw matrices, complete third-party database/network dumps, source proteomics tables, publication PDFs, or supplementary source files are redistributed; source-provider rights remain in force |
| Zenodo DOI | PENDING; no connection or deposit was attempted |

## Pre-commit inventory

The reviewed candidate contained the expected `README.md`, `scripts/`, `results/`, `environment/`, `PROVENANCE/`, `THIRD_PARTY_NOTICE.md`, `DATA_AND_DERIVED_OUTPUTS_NOTICE.md`, and `FILE_MANIFEST.csv`. The final scan covered author/contact data, telephone patterns, API keys, tokens, passwords, private keys, credentials, local absolute paths, prompts/conversations, Git internals, raw-expression/FASTQ formats, publication PDFs, and source supplementary files. No unresolved release blocker was found.

The 20 largest tracked files before the initial commit were:

| Rank | Relative path | Bytes |
|---:|---|---:|
| 1 | `results/04_phase4a/receptors/receiver_receptor_competence.csv` | 1,459,508 |
| 2 | `results/04_phase4a/integration/sender_receiver_evidence_matrix.csv` | 1,137,244 |
| 3 | `results/03_phase3/ligand_universe/frozen_ligand_universe.csv` | 794,612 |
| 4 | `results/03_phase3/baseline/baseline_sender_robustness.csv` | 593,931 |
| 5 | `results/03_phase3/sender/frozen_phase4_sender_candidates.csv` | 394,965 |
| 6 | `results/03_phase3/licensing/licensing_ligand_classification.csv` | 360,764 |
| 7 | `results/final_synthesis/final_sender_evidence.csv` | 122,413 |
| 8 | `results/04_phase4a/integration/phase4a_candidate_hierarchy.csv` | 84,330 |
| 9 | `results/04_phase4a/freeze/phase4_sender_scopes.csv` | 64,313 |
| 10 | `scripts/00_dataset_registry_phase0b/build_phase0b_audit.py` | 52,845 |
| 11 | `scripts/final_synthesis/build_final_synthesis.py` | 51,147 |
| 12 | `scripts/00_dataset_registry/build_audit.py` | 41,063 |
| 13 | `scripts/01_phase1a/build_phase1a_freeze.py` | 40,441 |
| 14 | `scripts/04_phase4b/build_phase4b_outputs.py` | 39,092 |
| 15 | `scripts/01_phase1b/run_phase1b_analysis.R` | 36,632 |
| 16 | `scripts/04_phase4a/run_phase4a_analysis.py` | 35,570 |
| 17 | `scripts/02_phase2a1/run_phase2a1_analysis.R` | 31,127 |
| 18 | `scripts/04_phase4b1/build_phase4b1_outputs.py` | 28,372 |
| 19 | `config/dataset_judgments.json` | 26,302 |
| 20 | `scripts/01_phase1a1/build_phase1a1_amendment.py` | 24,805 |

## GitHub authentication boundary

GitHub CLI was not installed, no `GH_TOKEN` or `GITHUB_TOKEN` was present, and no Git Credential Manager was available. Preparation therefore stops at a verified local Git repository. No password or token should be written into this repository.

After installing GitHub CLI and authenticating interactively, run from this repository:

```powershell
gh auth login
gh repo view hUCMSC-PE-paracrine-communication
```

If the second command confirms that no repository with that name exists in the intended account or organization, create a private repository and push:

```powershell
gh repo create hUCMSC-PE-paracrine-communication --private --description "Reproducible code and derived outputs for cross-dataset inference of hUC-MSC sender–receiver communication in preeclampsia." --source . --remote origin --push
git push origin v1.0.0
```

Inspect the private repository before any visibility change. Only after explicit human authorization should it be made public:

```powershell
gh repo edit hUCMSC-PE-paracrine-communication --visibility public
```

Only after separate human approval of release status should a GitHub Release be created:

```powershell
gh release create v1.0.0 --title "v1.0.0 — Manuscript reproducibility release" --notes "Reproducibility release accompanying the manuscript ‘Independent mapping of placental receiver states and hUC-MSC sender programs constrains putative paracrine communication in preeclampsia.’"
```

## Final gate

**LOCAL_PUBLIC_REPOSITORY_READY_FOR_GITHUB_AUTHENTICATION**

No scientific analysis was rerun, and no public GitHub repository or Zenodo record was created.
