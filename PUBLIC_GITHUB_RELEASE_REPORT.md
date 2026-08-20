# Public GitHub release preparation report

Preparation date: 2026-08-20

Repository name: `hUCMSC-PE-paracrine-communication`
Recommended description: Reproducible code and derived outputs for cross-dataset inference of hUC-MSC sender–receiver communication in preeclampsia.

## Repository status

| Field | Status |
|---|---|
| Local repository path | `release/public_repository_candidate` within the private project workspace; the user-specific absolute path is intentionally omitted from the public report |
| GitHub repository URL | `https://github.com/y15848102126-creator/hUCMSC-PE-paracrine-communication` |
| Visibility | **PUBLIC**; explicitly authorized after all private-repository checks passed |
| Current branch | `main`, tracking `origin/main` |
| Initial release commit SHA | `d12bfd3cb356ef5aebcd3cb469465e5eb68e8280` |
| Prepared release commit SHA | `834602d2049f301dbe10ecfef3f11908b90f6221` |
| Remote `main` immediately after the visibility-only operation | `240498efe74c5ee0bcb4ee779594c2cbebd6dc1a`; identical before and after the visibility change |
| README-correction base `main` SHA | `1ebb23cc05359d29ba21824facc297237824d4c0` |
| Existing tag status | Annotated `v1.0.0` remains unchanged: tag object `54e897952eec3afbdc060856b4541ef423d676c5`, peeled commit `834602d2049f301dbe10ecfef3f11908b90f6221` |
| Final archival tag status | Annotated `v1.0.1` unchanged after Release creation: tag object `087c96982ed5097ced7124c8ccfded87989c2ff7`, peeled commit `63bbded6d41e273e5880c2d745843b9c4e10b8ff` |
| GitHub Release status | CREATED and publicly accessible at `https://github.com/y15848102126-creator/hUCMSC-PE-paracrine-communication/releases/tag/v1.0.1` |
| GitHub Release title | `v1.0.1 — Manuscript reproducibility release` |
| GitHub Release tag/commit | Existing `v1.0.1`; peeled commit `63bbded6d41e273e5880c2d745843b9c4e10b8ff` |
| License | BSD 3-Clause for original project code only |
| Tracked files | 160, including this report and excluding `.git/` |
| Repository content size | 6,163,524 bytes, excluding `.git/` |
| Files excluded by `.gitignore` | 0 existing release-candidate files |
| Privacy/credential scan | PASS before Git initialization and unchanged during the push-only stage; no unresolved sensitive/private content detected |
| Manifest verification | PASS after regeneration; `FILE_MANIFEST.csv` excludes itself and covers every other tracked release file |
| Third-party redistribution status | No raw matrices, complete third-party database/network dumps, source proteomics tables, publication PDFs, or supplementary source files are redistributed; source-provider rights remain in force |
| Zenodo status | RELEASED; manuscript-facing version DOI `10.5281/zenodo.22026030` (concept DOI `10.5281/zenodo.22026029`) |
| Files changed by the visibility operation | NONE; the subsequent administrative update changed only this report and its `FILE_MANIFEST.csv` hash record |
| Files changed by GitHub Release creation | NONE; the subsequent administrative update changed only this report and its `FILE_MANIFEST.csv` hash record |

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

## GitHub push verification

- Git for Windows `2.53.0.windows.3` and GitHub CLI `2.97.0` were verified.
- GitHub CLI authentication was confirmed for `y15848102126-creator` without storing a token in the repository.
- An exact account-qualified lookup confirmed that `y15848102126-creator/hUCMSC-PE-paracrine-communication` did not exist before creation. No existing repository was overwritten.
- A new empty GitHub repository was created with `PRIVATE` visibility and the prespecified description.
- `origin` is exactly `https://github.com/y15848102126-creator/hUCMSC-PE-paracrine-communication.git`.
- The initial remote `main` SHA matched the prepared release commit `834602d2049f301dbe10ecfef3f11908b90f6221`.
- The remote annotated `v1.0.0` tag resolves to the same prepared release commit. The tag was not moved during this administrative report update.
- The initial remote history contained exactly the two independent public-repository commits, rooted at `d12bfd3cb356ef5aebcd3cb469465e5eb68e8280`; no parent private-project commit was present. The only subsequent `main` commit records this post-push report and its manifest hash.
- GitHub's rendered-README endpoint returned content containing the repository title and the authoritative receiver-method boundary.
- During the initial push stage, repository visibility was rechecked as `PRIVATE`; the default branch was `main`; the GitHub Release count was zero.
- The initial push stage performed no visibility change, GitHub Release creation, or Zenodo connection.

The SHA of the report-containing `main` commit is not embedded in this file because a commit cannot reproducibly contain its own SHA. It is verified by comparing local `HEAD`, `origin/main`, and the GitHub API after the report commit is pushed.

## Public-visibility verification

- Before the visibility change, the exact account/repository identity, private visibility, `main` default branch, sole `origin`, remote `main`, remote `v1.0.0`, independent root history, rendered README, required notices, zero GitHub Releases, manifest integrity, and privacy/credential scan all passed.
- The authenticated GitHub account changed only this repository's visibility from `PRIVATE` to `PUBLIC`; the repository was not deleted or recreated, and no force-push occurred.
- Immediately after the change, GitHub reported `PUBLIC` visibility and default branch `main`; unauthenticated repository access returned HTTP 200.
- Remote `main` remained `240498efe74c5ee0bcb4ee779594c2cbebd6dc1a` throughout the visibility-only operation.
- The annotated `v1.0.0` tag remained unchanged and peeled to `834602d2049f301dbe10ecfef3f11908b90f6221`.
- The README rendered successfully, and `LICENSE`, `THIRD_PARTY_NOTICE.md`, and `DATA_AND_DERIVED_OUTPUTS_NOTICE.md` remained accessible.
- GitHub Release count remained zero. No Zenodo interaction occurred.
- This post-verification report update is administrative and does not alter scientific outputs. Its commit advances `main` after the visibility checkpoint without moving `v1.0.0`.

## Administrative README correction and archival tag

- Human review identified one outdated publication-status sentence and the word `candidate` in the README title.
- The README title phrase was changed only from `public repository candidate v1.1` to `public repository v1.1`.
- The release-status wording was changed only to `PUBLIC REPOSITORY v1.1 — HUMAN RELEASE REVIEW PASSED`, followed by a concise notice that the archival release and Zenodo DOI remain pending.
- The only files changed in the correction commit were `README.md`, `PUBLIC_GITHUB_RELEASE_REPORT.md`, and the mechanically refreshed `FILE_MANIFEST.csv`.
- No file under `results/`, `scripts/`, `config/`, `environment/`, or `PROVENANCE/` changed.
- The normal commit message was `Update public repository release status`; no force-push occurred.
- Existing annotated tag `v1.0.0` was neither moved, deleted, recreated, nor force-updated.
- New annotated tag `v1.0.1`, with message `Public manuscript reproducibility release`, was created at the corrected public-release commit and pushed.
- The repository remained `PUBLIC`, the README rendered with the corrected status, and the GitHub Release count remained zero.
- No Zenodo interaction occurred.

## GitHub Release v1.0.1

- Immediately before Release creation, the repository was `PUBLIC`, GitHub Release count was zero, and existing annotated tag `v1.0.1` pointed to `63bbded6d41e273e5880c2d745843b9c4e10b8ff`.
- GitHub Release `v1.0.1 — Manuscript reproducibility release` was created with `--verify-tag`, so GitHub used the existing tag rather than creating or moving one.
- Release URL: `https://github.com/y15848102126-creator/hUCMSC-PE-paracrine-communication/releases/tag/v1.0.1`.
- The release is published, is neither a draft nor a prerelease, and contains no manually uploaded assets. GitHub's automatically generated source archives remain available.
- After Release creation, annotated tag `v1.0.1` retained object `087c96982ed5097ced7124c8ccfded87989c2ff7` and peeled commit `63bbded6d41e273e5880c2d745843b9c4e10b8ff`.
- Annotated tag `v1.0.0` also remained unchanged: object `54e897952eec3afbdc060856b4541ef423d676c5`, peeled commit `834602d2049f301dbe10ecfef3f11908b90f6221`.
- Repository visibility remained `PUBLIC`, and the Release URL returned HTTP 200 without authentication.
- GitHub Release creation changed no repository file. This later report/manifest update is administrative and does not modify scientific content or either release tag.
- Zenodo was not accessed directly. Automatic archival and DOI generation are pending through the user-enabled GitHub integration.

## Final gate

**PUBLIC_REPOSITORY_AND_ZENODO_ARCHIVE_RELEASED**

No scientific analysis was rerun. The repository and GitHub Release are public, and the immutable `v1.0.1` snapshot is archived at https://doi.org/10.5281/zenodo.22026030.

## Post-release DOI backfill

- The public README was updated administratively to display the released `v1.0.1` Zenodo archive.
- The manuscript-facing archival identifier is the version DOI `10.5281/zenodo.22026030`; the concept DOI was recorded only to prevent substitution.
- This post-release update advances `main` but does not move or recreate `v1.0.1`, create another GitHub Release, or create another Zenodo version.
- No scientific result, method, figure, threshold, candidate hierarchy, or evidence classification changed.
