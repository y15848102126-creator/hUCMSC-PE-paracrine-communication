#!/usr/bin/env python3
"""Create outcome-descriptive Phase 4B previews with Pillow only."""

from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/04_phase4b"
OUT = BASE / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FONT = ImageFont.load_default(size=18)
SMALL = ImageFont.load_default(size=14)
TITLE = ImageFont.load_default(size=24)


def canvas(width, height, title, subtitle=""):
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    draw.text((30, 22), title, font=TITLE, fill="#17202a")
    if subtitle:
        draw.text((30, 56), subtitle, font=SMALL, fill="#4d5656")
    return img, draw


def evidence_matrix():
    ev = pd.read_csv(BASE / "integration/phase4b_candidate_evidence_matrix.csv")
    cols = ["Topology", "Protein source", "Empirical reversal", "Mixed-direction risk"]
    img, draw = canvas(1160, 810, "Phase 4B independent evidence dimensions", "Categorical display only; no composite score")
    x0, y0, cw, rh = 300, 105, 190, 36
    for j, label in enumerate(cols):
        draw.text((x0 + j * cw + 8, y0 - 30), label, font=SMALL, fill="#17202a")
    protein_colors = {"NO_PROTEIN_SOURCE_EVIDENCE": "#eceff1", "HUCMSC_PROTEIN_CELL_ONLY": "#fff3cd", "OTHER_MSC_SECRETOME": "#ffe0b2", "HUCMSC_SECRETOME_DIRECT": "#80cbc4"}
    mixed_colors = {"SIGNED_EVIDENCE_UNRESOLVED": "#eceff1", "MIXED_DIRECTION_CONTEXT_DEPENDENT": "#ffcc80", "PREDOMINANTLY_DISEASE_CONCORDANT": "#ef9a9a"}
    for i, row in ev.iterrows():
        y = y0 + i * rh
        draw.text((30, y + 8), row.candidate, font=SMALL, fill="#17202a")
        plausible = row.PARACRINE_TOPOLOGY not in {"MEMBRANE_ASSOCIATED", "INTRACELLULAR_OR_QUESTIONABLE_PARACRINE", "UNCERTAIN"}
        reversal = int(str(row.EMPIRICAL_SIGNED_PERTURBATION).split(";")[0].split("=")[1]) > 0
        cells = ["#90caf9" if plausible else "#ef9a9a", protein_colors.get(row.HUCMSC_PROTEIN_SOURCE, "#eceff1"), "#80cbc4" if reversal else "#eceff1", mixed_colors.get(row.MIXED_DIRECTION_RISK, "#eceff1")]
        labels = ["plausible" if plausible else "weak", row.HUCMSC_PROTEIN_SOURCE.replace("_", " ")[:18], "yes" if reversal else "none", row.MIXED_DIRECTION_RISK.replace("_", " ")[:18]]
        for j, (color, label) in enumerate(zip(cells, labels)):
            x = x0 + j * cw
            draw.rectangle((x, y, x + cw - 4, y + rh - 4), fill=color, outline="#ffffff")
            draw.text((x + 7, y + 8), label, font=SMALL, fill="#17202a")
    img.save(OUT / "phase4b_evidence_dimension_matrix.png")


def bar_chart(name, title, subtitle, items, colors, x_max):
    img, draw = canvas(1200, 500, title, subtitle)
    x0, y0, width, rh = 470, 110, 650, 70
    for i, (label, value) in enumerate(items):
        y = y0 + i * rh
        draw.text((30, y + 14), label.replace("_", " "), font=SMALL, fill="#17202a")
        bar = int(width * value / max(x_max, 1))
        draw.rectangle((x0, y, x0 + bar, y + 38), fill=colors[i], outline=colors[i])
        draw.text((x0 + bar + 10, y + 10), str(value), font=FONT, fill="#17202a")
    img.save(OUT / name)


def perturbation_counts():
    q = pd.read_csv(BASE / "perturbation/empirical_signed_perturbation_evidence.csv")
    order = ["EMPIRICAL_REVERSAL_SUPPORTED", "EMPIRICAL_DISEASE_CONCORDANT", "EMPIRICAL_SIGN_CONFLICT", "NO_EMPIRICAL_SIGNED_EVIDENCE"]
    counts = q.empirical_signed_classification.value_counts()
    items = [(x, int(counts.get(x, 0))) for x in order]
    bar_chart("phase4b_empirical_perturbation_summary.png", "Independent module-matched perturbation evidence", "Frozen Phase 4A Tier A ligand-module axes (n=38)", items, ["#2a9d8f", "#e76f51", "#f4a261", "#9aa0a6"], 38)


def candidate_states():
    c = pd.read_csv(BASE / "integration/phase4b_candidate_classification.csv")
    labels = ["TRIANGULATED_HIGH_PRIORITY", "TRIANGULATED_CONTEXT_DEPENDENT", "PERTURBATION_SUPPORTED_BUT_SOURCE_UNCONFIRMED", "PROTEIN_SUPPORTED_BUT_DIRECTION_UNRESOLVED", "PARTIAL_EXTERNAL_EVIDENCE", "COMPUTATIONAL_ONLY", "BIOPHYSICALLY_WEAK_PARACRINE"]
    primary = c.primary_classification.value_counts()
    values = [int(primary.get(x, 0)) if x not in {"COMPUTATIONAL_ONLY", "BIOPHYSICALLY_WEAK_PARACRINE"} else int((c[x] == "YES").sum()) for x in labels]
    img, draw = canvas(1260, 720, "Deterministic Phase 4B classification outcome", "COMPUTATIONAL_ONLY and BIOPHYSICALLY_WEAK are non-exclusive flags")
    x0, y0, width, rh = 550, 105, 610, 80
    colors = ["#006d77", "#83c5be", "#457b9d", "#8ecae6", "#adb5bd", "#b0bec5", "#d1495b"]
    for i, (label, value) in enumerate(zip(labels, values)):
        y = y0 + i * rh
        draw.text((30, y + 12), label.replace("_", " "), font=SMALL, fill="#17202a")
        bar = int(width * value / 17)
        draw.rectangle((x0, y, x0 + bar, y + 42), fill=colors[i])
        draw.text((x0 + bar + 10, y + 10), str(value), font=FONT, fill="#17202a")
    img.save(OUT / "phase4b_candidate_classification.png")


if __name__ == "__main__":
    evidence_matrix()
    perturbation_counts()
    candidate_states()
    print("PHASE4B_FIGURES_OK: 3 previews")
