#!/usr/bin/env python3
"""Attach explicit source fields to every derived Phase 2A analytical CSV."""
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/"results/02_phase2a"
DATA="https://doi.org/10.6084/m9.figshare.23264102.v1"
PAPER="https://doi.org/10.1016/j.medj.2023.07.005"
MSIG="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/"
COLLECTRI="https://github.com/saezlab/CollecTRI|https://omnipathdb.org/"

def amend(path: Path, url: str, accession: str) -> None:
    frame=pd.read_csv(path,low_memory=False)
    frame["source_url"]=url; frame["source_accession"]=accession
    frame.to_csv(path,index=False,lineterminator="\n")

for path in (OUT/"DE").glob("*.csv"): amend(path,f"{DATA}|{PAPER}","Figshare:23264102.v1:file41003240")
for path in (OUT/"programs").glob("*.csv"): amend(path,f"{DATA}|{MSIG}","Figshare:23264102.v1;MSigDB:2026.1.Hs")
amend(OUT/"regulons/cellstate_regulon_activity.csv",f"{DATA}|{COLLECTRI}","Figshare:23264102.v1;CollecTRI:OmniPath_snapshot_2026-08-09")
amend(OUT/"qc/pseudobulk_mds_coordinates.csv",DATA,"Figshare:23264102.v1:file41003240")
amend(OUT/"qc/phase2a_qc_summary.csv",f"{DATA}|config/phase2a_analysis.json","Figshare:23264102.v1;PHASE2A_CONFIG")
print("Phase 2A provenance fields added")
