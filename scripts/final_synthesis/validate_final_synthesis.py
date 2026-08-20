#!/usr/bin/env python3
"""Validate frozen final-synthesis invariants without rerunning analyses."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EXPECTED_ROWS = {
    "phase_registry.csv": 13,
    "final_receiver_evidence.csv": 11,
    "final_sender_evidence.csv": 214,
    "final_17_candidate_evidence.csv": 17,
    "negative_results_register.csv": 8,
    "claim_ceiling_matrix.csv": 12,
    "figure_source_registry.csv": 16,
    "supplementary_table_registry.csv": 16,
    "experimental_validation_roadmap.csv": 8,
    "final_synthesis_risk_flags.csv": 12,
}


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    out = root / "results" / "final_synthesis"

    tables = {}
    for name, n in EXPECTED_ROWS.items():
        path = out / name
        assert path.exists(), path
        df = read(path)
        assert len(df) == n, (name, len(df), n)
        assert "source_url" in df.columns and df["source_url"].str.len().gt(0).all(), name
        tables[name] = df

    phases = tables["phase_registry.csv"].set_index("phase")
    assert phases.loc["1B", "result_status"] == "FROZEN_NEGATIVE"
    assert phases.loc["2A", "result_status"] == "LEGACY_COUNT_MODEL_DISCOVERY"
    assert phases.loc["4B", "result_status"] == "FROZEN_SUPERSEDED_PROTEIN_DIMENSION"

    receiver = tables["final_receiver_evidence.csv"].set_index("module")
    assert receiver.loc["PROGRAM_MODULE_07", "final_receiver_level"] == "R1"
    assert receiver.loc["PROGRAM_MODULE_08", "final_receiver_level"] == "HOLD"
    assert receiver["bulk_robust_evidence"].eq("NO").all()

    sender = tables["final_sender_evidence.csv"]
    assert sender["sender_evidence_level"].value_counts().to_dict() == {"S2": 176, "S1": 38}
    assert sender["P1_PARACRINE_CORE"].eq("YES").sum() == 148
    assert sender["P2_EXTRACELLULAR_EXTENDED"].eq("YES").sum() == 190
    assert sender["P3_FULL_LR_SENSITIVITY"].eq("YES").sum() == 214
    assert sender["huc_wj_msc_ev_protein_evidence"].eq("NOT_SYSTEMATICALLY_AUDITED").sum() == 197

    candidates = tables["final_17_candidate_evidence.csv"].set_index("candidate")
    assert candidates["final_deterministic_category"].eq("TRIANGULATED_HIGH_PRIORITY").sum() == 1
    enpp1 = candidates.loc["ENPP1"]
    assert enpp1["final_frozen_status"] == "MECHANICALLY_TRIANGULATED_LEAD_HYPOTHESIS"
    assert enpp1["context_dependent"] == "YES"
    assert enpp1["ev_only_protein_source"] == "YES"
    assert enpp1["non_placental_empirical_reversal"] == "YES"
    assert enpp1["direct_pe_experimental_validation"] == "NO"

    neg = tables["negative_results_register.csv"]
    assert {"NEG01", "NEG02", "NEG03", "NEG04", "NEG05", "NEG06", "NEG07", "NEG08"} == set(neg["negative_result_id"])

    claims = tables["claim_ceiling_matrix.csv"]
    assert not claims["maximum_allowed_wording"].str.contains("validated therapeutic factor", case=False).any()
    assert not claims["supported_claim"].str.contains("hUC-MSCs reverse", case=False).any()

    fig = tables["figure_source_registry.csv"]
    assert all((root / rel).exists() for rel in fig["source_file"])
    supp = tables["supplementary_table_registry.csv"]
    assert not supp["source_row_counts"].str.contains("MISSING").any()
    for sources in supp["source_files"]:
        assert all((root / rel).exists() for rel in sources.split("|")), sources

    docs = [
        "FINAL_SYNTHESIS_REPORT.md", "MANUSCRIPT_RESULTS_BLUEPRINT.md",
        "MANUSCRIPT_DISCUSSION_BLUEPRINT.md", "FIGURE_PLAN_FINAL.md", "SUPPLEMENTARY_PLAN_FINAL.md",
    ]
    for name in docs:
        text = (root / "docs" / name).read_text(encoding="utf-8")
        assert len(text) > 1000, name
    report = (root / "docs" / "FINAL_SYNTHESIS_REPORT.md").read_text(encoding="utf-8")
    assert "READY_WITH_RESTRICTIONS" in report
    assert "MECHANICALLY_TRIANGULATED_LEAD_HYPOTHESIS" in report
    assert "validated therapeutic factor" in report and "Prohibited wording" in report

    print("PASS: final synthesis tables, frozen invariants, source paths, and documents validated")


if __name__ == "__main__":
    main()
