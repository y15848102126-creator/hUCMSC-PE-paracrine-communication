# Third-party resource retrieval

## MSigDB 2026.1.Hs

Retrieve Hallmark, Reactome, and GO Biological Process symbol GMT files from the official MSigDB release endpoint recorded in `scripts/02_phase2a/download_phase2a_resources.py`. Provider access and license terms apply. Reconstruct frozen hypotheses with `scripts/02_phase2b/freeze_phase2b_hypotheses.py` and verify membership sizes and SHA-256 values against the sanitized hypothesis tables.

## NicheNet and OmniPath/CollecTRI

Use `scripts/04_phase4a/prepare_nichenet_resources.R`, `scripts/04_phase4a/download_phase4a_resources.py`, and `scripts/02_phase2a/download_phase2a_resources.py` to retrieve resources from the original providers. Exact snapshot identifiers and hashes are retained in configurations and derived registries; no complete network dump is redistributed.

## Public transcriptomic and proteomic data

Use persistent accessions in `submission/RBE_Final_Dataset_Reference_Registry.csv`. GEO, Figshare, GitHub, and ProteomeXchange/PRIDE source files are not repackaged.
