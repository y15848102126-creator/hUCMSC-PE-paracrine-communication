# Zenodo release instructions

No Zenodo connection or deposit has been created by this repository-preparation step. Perform the following only after the GitHub repository has passed human inspection and its public visibility and release status have been explicitly approved.

1. Sign in to [Zenodo](https://zenodo.org/) with the account that will own the record.
2. Open the GitHub integration from the Zenodo account settings and authorize the appropriate GitHub organization or account.
3. Locate `hUCMSC-PE-paracrine-communication` in the Zenodo GitHub repository list and enable archiving for that repository.
4. On GitHub, verify that the repository content, metadata, license boundary, third-party notices, and visibility are final.
5. Push or verify the annotated `v1.0.0` tag, then create the GitHub Release titled `v1.0.0 — Manuscript reproducibility release` using that exact tag.
6. Use this release description: `Reproducibility release accompanying the manuscript “Independent mapping of placental receiver states and hUC-MSC sender programs constrains putative paracrine communication in preeclampsia.”`
7. Wait for Zenodo to archive GitHub release `v1.0.0`; do not create a competing manual upload for the same release.
8. Review the draft Zenodo record before publication. Verify title, author names and order, affiliations, ORCID identifiers, description, keywords, version, publication date, license scope, related identifiers, and repository URL. Do not imply that BSD-3-Clause covers third-party data or database resources.
9. Publish the Zenodo record only after author approval. Record the version DOI assigned to `v1.0.0` and, if Zenodo supplies one, the concept DOI covering all versions.
10. Add the approved DOI metadata to the manuscript data/code-availability statement and repository README in a later, separately reviewed update.

Before publication, confirm that the Zenodo archive contains no raw matrices, controlled data, credentials, private notes, source publication PDFs or supplements, or third-party resources that must be retrieved from their original providers.
