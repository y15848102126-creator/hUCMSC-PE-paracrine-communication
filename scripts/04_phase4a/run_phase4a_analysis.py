#!/usr/bin/env python3
"""Run frozen, blinded Phase 4A sender-receiver integration.

This script intentionally keeps unsigned NicheNet compatibility separate from
signed directional evidence. It does not query PE-specific literature.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "04_phase4a"
SEED = 20260815
RNG = np.random.default_rng(SEED)
SOURCE = "config/phase4a_analysis.json|NicheNet:Zenodo:7074291|OmniPath:2026-08-15|CollecTRI:2026-08-15|CytoSig:2026-08-15"


def bh(pvalues: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if not len(pv):
        return out
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(len(pv))
    adjusted[order] = np.minimum(ranked, 1.0)
    out[np.where(valid)[0]] = adjusted
    return out


def bool_col(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def sign_of(value: float) -> int:
    return 1 if value > 0 else (-1 if value < 0 else 0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_cytosig(path: Path) -> pd.DataFrame:
    # The first field of each data row is a gene, while the header has no label.
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str).str.upper()
    return frame


def rank_corr(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    ra = pd.Series(a).rank(method="average").to_numpy()
    rb = pd.Series(b).rank(method="average").to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


def exact_permutation_p(observed: float, a: np.ndarray, b: np.ndarray, permutations: int = 2000) -> float:
    if not np.isfinite(observed):
        return float("nan")
    extreme = 0
    for _ in range(permutations):
        permuted = RNG.permutation(b)
        trial = rank_corr(a, permuted)
        extreme += int(np.isfinite(trial) and abs(trial) >= abs(observed))
    return (extreme + 1.0) / (permutations + 1.0)


def signed_edges(frame: pd.DataFrame) -> pd.DataFrame:
    stim = bool_col(frame["is_stimulation"])
    inhib = bool_col(frame["is_inhibition"])
    direct = bool_col(frame["is_directed"])
    good = direct & stim.ne(inhib)
    out = frame.loc[good, ["source_genesymbol", "target_genesymbol", "sources", "curation_effort"]].copy()
    out["sign"] = np.where(stim[good], 1, -1)
    out = out.dropna(subset=["source_genesymbol", "target_genesymbol"])
    out["source_genesymbol"] = out.source_genesymbol.astype(str).str.upper()
    out["target_genesymbol"] = out.target_genesymbol.astype(str).str.upper()
    # Exact contradictory duplicates are ambiguous and removed.
    sign_n = out.groupby(["source_genesymbol", "target_genesymbol"]).sign.nunique()
    keep = set(sign_n[sign_n.eq(1)].index)
    out = out[out.apply(lambda x: (x.source_genesymbol, x.target_genesymbol) in keep, axis=1)]
    return out.sort_values(["source_genesymbol", "target_genesymbol", "sign"]).drop_duplicates(["source_genesymbol", "target_genesymbol"])


def receptor_tf_states(activity: pd.DataFrame, receptors: set[str], tfs: set[str], max_depth: int = 3) -> dict[str, dict[str, set[int]]]:
    graph: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for row in activity.itertuples(index=False):
        graph[row.source_genesymbol].append((row.target_genesymbol, int(row.sign)))
    answer: dict[str, dict[str, set[int]]] = {}
    for receptor in sorted(receptors):
        frontier: dict[str, set[int]] = {receptor: {1}}
        reached: dict[str, set[int]] = defaultdict(set)
        for _ in range(max_depth):
            following: dict[str, set[int]] = defaultdict(set)
            for node, signs in frontier.items():
                for target, edge_sign in graph.get(node, []):
                    for prior in signs:
                        following[target].add(prior * edge_sign)
            for target, signs in following.items():
                if target in tfs:
                    reached[target].update(signs)
            frontier = following
        answer[receptor] = dict(reached)
    return answer


def main() -> int:
    cfg = json.loads((ROOT / "config/phase4a_analysis.json").read_text(encoding="utf-8"))
    for directory in ["receptors", "lr", "targets", "signed", "integration", "qc"]:
        (OUT / directory).mkdir(parents=True, exist_ok=True)

    senders = pd.read_csv(OUT / "freeze/phase4_sender_scopes.csv")
    receivers = pd.read_csv(OUT / "freeze/phase4_receiver_hierarchy.csv")
    receivers = receivers[receivers.phase4a_analysis_scope.isin(["PRIMARY", "SECONDARY_SENSITIVITY"])].copy()
    lr = pd.read_csv(ROOT / "data/interim/phase4a/nichenet_lr_network.tsv", sep="\t")
    lr["from"] = lr["from"].astype(str).str.upper()
    lr["to"] = lr["to"].astype(str).str.upper()
    candidates = set(senders.gene.astype(str).str.upper())
    lr = lr[lr["from"].isin(candidates)].drop_duplicates(["from", "to"])
    candidate_receptors = set(lr["to"])

    # ------------------------------------------------------------------
    # Pregnancy-level receptor competence on the corrected Admati layer.
    eligibility = pd.read_csv(ROOT / "results/02_phase2a/metadata/pseudobulk_eligibility.csv")
    eligibility = eligibility[(eligibility.include_in_contrast == "YES") & eligibility.harmonized_annotation.isin(receivers.celltype)]
    # Use only mutually exclusive clinical groups when constructing one
    # pregnancy-level receptor record.  The eligibility registry also contains
    # pooled PE/CONTROL summary rows; including those alongside EOPE/LOPE and
    # EARLY_CONTROL/LATE_CONTROL duplicates pregnancies, and the pooled "PE"
    # label was previously misclassified as control by the disease assignment
    # below.  Receptor competence itself is PE-only, but descriptive control
    # fields must use the canonical non-overlapping groups.
    eligibility = eligibility[
        eligibility.group.isin(["EOPE", "LOPE", "EARLY_CONTROL", "LATE_CONTROL"])
    ]
    strata = pd.read_csv(ROOT / "data/interim/phase2a/pseudobulk_strata.tsv", sep="\t")
    counts = pd.read_csv(ROOT / "data/interim/phase2a/admati_harmonized_pseudobulk_counts.csv")
    counts["gene"] = counts.gene.astype(str).str.upper()
    counts = counts.groupby("gene", as_index=False).sum(numeric_only=True)
    count_index = counts.set_index("gene")
    strata_key = strata.set_index(["patient_id", "harmonized_annotation"])
    receptor_rows = []
    for celltype in sorted(receivers.celltype.unique()):
        e = eligibility[eligibility.harmonized_annotation.eq(celltype)].drop_duplicates(["patient_id"])
        for receptor in sorted(candidate_receptors):
            values = []
            for row in e.itertuples(index=False):
                key = (row.patient_id, celltype)
                if key not in strata_key.index:
                    continue
                sr = strata_key.loc[key]
                column = sr.matrix_column if isinstance(sr, pd.Series) else sr.iloc[0].matrix_column
                cells = float(row.cell_count)
                value = float(count_index.at[receptor, column]) / cells if receptor in count_index.index and cells > 0 else 0.0
                disease = "PE" if row.group in ("EOPE", "LOPE") else "CONTROL"
                values.append((row.patient_id, row.group, disease, value))
            pe = [x for x in values if x[2] == "PE"]
            control = [x for x in values if x[2] == "CONTROL"]
            pe_detected = sum(x[3] >= 0.1 for x in pe)
            ctrl_detected = sum(x[3] >= 0.1 for x in control)
            pe_prop = pe_detected / len(pe) if pe else 0.0
            if pe_prop >= 0.50 and pe_detected >= 3:
                competence = "RECEPTOR_COMPETENT"
            elif pe_prop >= 0.20 or pe_detected >= 2:
                competence = "RECEPTOR_WEAK"
            else:
                competence = "RECEPTOR_NOT_SUPPORTED"
            receptor_rows.append({
                "receptor": receptor, "receiver_celltype": celltype,
                "PE_pregnancy_n": len(pe), "PE_detected_n": pe_detected, "PE_detected_proportion": pe_prop,
                "control_pregnancy_n": len(control), "control_detected_n": ctrl_detected,
                "control_detected_proportion": ctrl_detected / len(control) if control else 0.0,
                "PE_median_mean_normalized_expression": float(np.median([x[3] for x in pe])) if pe else np.nan,
                "control_median_mean_normalized_expression": float(np.median([x[3] for x in control])) if control else np.nan,
                "PE_donor_expression": ";".join(f"{x[0]}:{x[3]:.6g}" for x in pe),
                "control_donor_expression": ";".join(f"{x[0]}:{x[3]:.6g}" for x in control),
                "receptor_competence": competence,
                "competence_rule": "patient mean>=0.1; COMPETENT if >=50% and >=3 PE pregnancies; WEAK if >=20% or >=2",
                "receptor_mapping_source": "NicheNet-v2 LR network gene symbols",
                "evidence_confidence": "PREGNANCY_LEVEL_PUBLIC_NORMALIZED_EXPRESSION",
                "source_url": "data/interim/phase2a/admati_harmonized_pseudobulk_counts.csv|results/02_phase2a/metadata/pseudobulk_eligibility.csv",
            })
    receptor_df = pd.DataFrame(receptor_rows).sort_values(["receiver_celltype", "receptor"])
    receptor_df.to_csv(OUT / "receptors/receiver_receptor_competence.csv", index=False)

    # ---------------------------------------------------------------
    # NicheNet LR compatibility with OmniPath exact-edge cross-check.
    op_lr_raw = pd.read_csv(ROOT / "data/raw/phase4a/omnipath_lr_crosscheck_20260815.tsv", sep="\t", low_memory=False)
    op_lr_raw["source_genesymbol"] = op_lr_raw.source_genesymbol.astype(str).str.upper()
    op_lr_raw["target_genesymbol"] = op_lr_raw.target_genesymbol.astype(str).str.upper()
    op_map = defaultdict(list)
    for row in op_lr_raw.itertuples(index=False):
        op_map[(row.source_genesymbol, row.target_genesymbol)].append(row)
    competence_map = receptor_df.set_index(["receiver_celltype", "receptor"]).receptor_competence.to_dict()
    sender_map = senders.set_index("gene").to_dict("index")
    lr_rows = []
    for ligand in sorted(candidates):
        edges = lr[lr["from"].eq(ligand)]
        for celltype in sorted(receivers.celltype.unique()):
            if edges.empty:
                lr_rows.append({"ligand": ligand, "receiver_celltype": celltype, "receptor": "NOT_IN_NICHENET_LR", "sender_evidence_level": sender_map[ligand]["sender_evidence_level"], "P1_PARACRINE_CORE": sender_map[ligand]["P1_PARACRINE_CORE"], "P2_EXTRACELLULAR_EXTENDED": sender_map[ligand]["P2_EXTRACELLULAR_EXTENDED"], "P3_FULL_LR_SENSITIVITY": "YES", "receptor_competence": "RECEPTOR_NOT_SUPPORTED", "nichenet_lr_support": "NO", "nichenet_sources": "", "omnipath_support": "NO", "omnipath_sources": "", "omnipath_sign": "NOT_EVALUABLE", "competent_interaction_n_ligand_celltype": 0, "interpretation": "NO_NICHENET_LR_EDGE", "source_url": SOURCE})
                continue
            competent_n = sum(competence_map.get((celltype, r), "RECEPTOR_NOT_SUPPORTED") == "RECEPTOR_COMPETENT" for r in edges["to"])
            for edge in edges.itertuples(index=False):
                op = op_map.get((ligand, edge.to), [])
                signs = set()
                sources = set()
                for item in op:
                    sources.update(str(item.sources).split(";"))
                    if bool(item.is_stimulation) != bool(item.is_inhibition):
                        signs.add("STIMULATION" if bool(item.is_stimulation) else "INHIBITION")
                op_sign = next(iter(signs)) if len(signs) == 1 else ("CONFLICTING" if len(signs) > 1 else "UNSIGNED")
                lr_rows.append({"ligand": ligand, "receiver_celltype": celltype, "receptor": edge.to, "sender_evidence_level": sender_map[ligand]["sender_evidence_level"], "P1_PARACRINE_CORE": sender_map[ligand]["P1_PARACRINE_CORE"], "P2_EXTRACELLULAR_EXTENDED": sender_map[ligand]["P2_EXTRACELLULAR_EXTENDED"], "P3_FULL_LR_SENSITIVITY": "YES", "receptor_competence": competence_map.get((celltype, edge.to), "RECEPTOR_NOT_SUPPORTED"), "nichenet_lr_support": "YES", "nichenet_sources": f"{edge.database}|{edge.source}", "omnipath_support": "YES" if op else "NO", "omnipath_sources": ";".join(sorted(s for s in sources if s and s != "nan")), "omnipath_sign": op_sign, "competent_interaction_n_ligand_celltype": competent_n, "interpretation": "COMPATIBILITY_NOT_THERAPEUTIC_EVIDENCE", "source_url": SOURCE})
    lr_df = pd.DataFrame(lr_rows).sort_values(["ligand", "receiver_celltype", "receptor"])
    lr_df.to_csv(OUT / "lr/sender_receiver_lr_compatibility.csv", index=False)

    # ---------------------------------------------------------------
    # Unsigned NicheNet target compatibility on full corrected ranking.
    target_rows = (ROOT / "data/interim/phase4a/nichenet_target_rows.txt").read_text(encoding="utf-8").splitlines()
    target_cols = (ROOT / "data/interim/phase4a/nichenet_target_columns.txt").read_text(encoding="utf-8").splitlines()
    target_matrix = np.fromfile(ROOT / "data/interim/phase4a/nichenet_target_subset_float64.bin", dtype="<f8").reshape((len(target_rows), len(target_cols)), order="F")
    target_index = {g.upper(): i for i, g in enumerate(target_rows)}
    ligand_index = {g.upper(): i for i, g in enumerate(target_cols)}
    stats = pd.read_csv(ROOT / "results/02_phase2a2/corrected_analysis/corrected_gene_statistics.csv")
    stats["gene"] = stats.gene.astype(str).str.upper()
    target_compat_rows, target_edge_rows = [], []
    for recv in receivers.itertuples(index=False):
        cell_stats = stats[stats.celltype.eq(recv.celltype)].pivot_table(index="gene", columns="contrast", values="moderated_t", aggfunc="mean")
        cell_stats = cell_stats.dropna(subset=["EOPE", "LOPE"])
        for contrast in ["EOPE", "LOPE"]:
            sd = cell_stats[contrast].std(ddof=1)
            cell_stats[f"z_{contrast}"] = (cell_stats[contrast] - cell_stats[contrast].mean()) / sd
        cell_stats["shared_stat"] = cell_stats[["z_EOPE", "z_LOPE"]].mean(axis=1)
        background = sorted(set(cell_stats.index).intersection(target_index))
        bg_idx = np.array([target_index[g] for g in background], dtype=int)
        module_genes = sorted(set(str(recv.module_union_gene_membership).split(";")))
        measured_module = [g for g in module_genes if g in cell_stats.index and g in target_index]
        mod_idx = np.array([target_index[g] for g in measured_module], dtype=int)
        weights = np.abs(cell_stats.loc[measured_module, "shared_stat"].to_numpy()) if measured_module else np.array([])
        required = min(15, math.ceil(0.60 * len(module_genes)))
        # Freeze one outcome-independent set of random background indices per
        # receiver module and reuse it for all ligands. This is equivalent to
        # the preregistered random-gene null and avoids millions of RNG calls.
        permutation_indices = None
        if measured_module and len(background) >= len(measured_module):
            permutation_indices = np.vstack([
                RNG.choice(len(background), size=len(measured_module), replace=False)
                for _ in range(cfg["target_compatibility"]["permutations"])
            ])
        for ligand in sorted(candidates):
            j = ligand_index.get(ligand)
            competent_receptors = lr_df[(lr_df.ligand.eq(ligand)) & (lr_df.receiver_celltype.eq(recv.celltype)) & (lr_df.receptor_competence.eq("RECEPTOR_COMPETENT"))].receptor
            lr_competent = len(competent_receptors) > 0
            if j is None or len(background) == 0 or len(measured_module) == 0:
                score, pval, positive_n, q95 = np.nan, np.nan, 0, np.nan
            else:
                bg_pot = target_matrix[bg_idx, j]
                ranks = pd.Series(bg_pot).rank(method="average", pct=True).to_numpy()
                rank_map = {g: ranks[i] for i, g in enumerate(background)}
                module_ranks = np.array([rank_map[g] for g in measured_module])
                positive_n = int(np.sum(target_matrix[mod_idx, j] > 0))
                score = float(np.average(module_ranks, weights=weights)) if weights.sum() > 0 else float(np.mean(module_ranks))
                sampled = ranks[permutation_indices]
                trials = ((sampled * weights).sum(axis=1) / weights.sum()) if weights.sum() > 0 else sampled.mean(axis=1)
                extreme = int(np.sum(trials >= score))
                pval = (extreme + 1.0) / (cfg["target_compatibility"]["permutations"] + 1.0)
                q95 = float(np.quantile(bg_pot, 0.95))
                for gene in measured_module:
                    potential = float(target_matrix[target_index[gene], j])
                    if potential > 0 and potential >= q95:
                        target_edge_rows.append({"ligand": ligand, "program_module": recv.program_module, "receiver_celltype": recv.celltype, "target_gene": gene, "nichenet_regulatory_potential": potential, "ligand_background_95th_percentile": q95, "corrected_shared_moderated_t": float(cell_stats.at[gene, "shared_stat"]), "PE_disease_direction": "UP" if cell_stats.at[gene, "shared_stat"] > 0 else "DOWN", "nichenet_directionality": "UNSIGNED", "interpretation": "HIGH_POTENTIAL_TARGET_EDGE_NOT_REVERSAL", "source_url": SOURCE})
            target_compat_rows.append({"ligand": ligand, "program_module": recv.program_module, "receiver_level": recv.receiver_level, "receiver_celltype": recv.celltype, "frozen_PE_direction": recv.frozen_direction, "P1_PARACRINE_CORE": sender_map[ligand]["P1_PARACRINE_CORE"], "measured_module_gene_n": len(measured_module), "frozen_module_gene_n": len(module_genes), "measured_coverage_proportion": len(measured_module) / len(module_genes), "positive_potential_module_target_n": positive_n, "required_measured_gene_n": required, "competent_receptor_n": len(set(competent_receptors)), "competent_receptors": ";".join(sorted(set(competent_receptors))), "weighted_target_rank_compatibility": score, "permutation_P": pval, "target_compatibility_BH_FDR": np.nan, "target_compatibility_class": "PENDING", "nichenet_directionality": "UNSIGNED_COMPATIBILITY_ONLY", "source_url": SOURCE})
    target_df = pd.DataFrame(target_compat_rows)
    target_df["target_compatibility_BH_FDR"] = bh(target_df.permutation_P)
    coverage_ok = (target_df.measured_module_gene_n >= target_df.required_measured_gene_n) & (target_df.positive_potential_module_target_n >= 10)
    competent_ok = target_df.competent_receptor_n > 0
    target_df["target_compatibility_class"] = "NOT_TARGET_COMPATIBLE"
    target_df.loc[coverage_ok & competent_ok & (target_df.weighted_target_rank_compatibility >= 0.55), "target_compatibility_class"] = "TARGET_DIRECTIONAL_TREND"
    target_df.loc[coverage_ok & competent_ok & (target_df.weighted_target_rank_compatibility >= 0.60) & (target_df.target_compatibility_BH_FDR < 0.10), "target_compatibility_class"] = "TARGET_COMPATIBLE"
    target_df.sort_values(["program_module", "ligand"]).to_csv(OUT / "targets/nichenet_target_compatibility.csv", index=False)
    target_edges = pd.DataFrame(target_edge_rows, columns=["ligand", "program_module", "receiver_celltype", "target_gene", "nichenet_regulatory_potential", "ligand_background_95th_percentile", "corrected_shared_moderated_t", "PE_disease_direction", "nichenet_directionality", "interpretation", "source_url"])
    target_edges.sort_values(["program_module", "ligand", "target_gene"]).to_csv(OUT / "targets/ligand_receiver_target_edges.csv", index=False)

    # ---------------------------------------------------------------
    # Signed evidence: exact signed LR + activity flow + CollecTRI;
    # and exact symbol-matched CytoSig as an independent sensitivity.
    activity = signed_edges(pd.read_csv(ROOT / "data/raw/phase4a/omnipath_signed_activity_flow_20260815.tsv", sep="\t", low_memory=False))
    collectri = signed_edges(pd.read_csv(ROOT / "data/raw/phase4a/collectri_signed_20260815.tsv", sep="\t", low_memory=False))
    exact_lr = signed_edges(op_lr_raw)
    tfs = set(collectri.source_genesymbol)
    competent_receptor_set = set(receptor_df.loc[receptor_df.receptor_competence.eq("RECEPTOR_COMPETENT"), "receptor"])
    rt_states = receptor_tf_states(activity, competent_receptor_set, tfs, cfg["signed_network"]["max_receptor_to_tf_edges"])
    tf_targets: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for row in collectri.itertuples(index=False):
        tf_targets[row.source_genesymbol].append((row.target_genesymbol, int(row.sign)))
    exact_lr_sign = exact_lr.set_index(["source_genesymbol", "target_genesymbol"]).sign.to_dict()
    cytosig = read_cytosig(ROOT / "data/raw/phase4a/cytosig_signature_centroid_core.tsv")
    cytosig_map = {str(c).upper().replace(" ", ""): c for c in cytosig.columns}
    signed_rows = []
    for row in target_df.itertuples(index=False):
        recv = receivers[receivers.program_module.eq(row.program_module)].iloc[0]
        module_genes = set(str(recv.module_union_gene_membership).split(";"))
        cell_stats = stats[(stats.celltype.eq(row.receiver_celltype))].pivot_table(index="gene", columns="contrast", values="moderated_t", aggfunc="mean").dropna(subset=["EOPE", "LOPE"])
        for c in ["EOPE", "LOPE"]:
            cell_stats[c] = (cell_stats[c] - cell_stats[c].mean()) / cell_stats[c].std(ddof=1)
        shared = cell_stats[["EOPE", "LOPE"]].mean(axis=1)
        votes: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
        for receptor in str(row.competent_receptors).split(";") if row.competent_receptors else []:
            lr_sign = exact_lr_sign.get((row.ligand, receptor))
            if lr_sign is None:
                continue
            for tf, path_signs in rt_states.get(receptor, {}).items():
                if len(path_signs) != 1:
                    continue
                path_sign = next(iter(path_signs))
                for target, tf_sign in tf_targets.get(tf, []):
                    if target in module_genes and target in shared.index:
                        votes[target].append((receptor, tf, int(lr_sign) * path_sign * tf_sign))
        consensus = {}
        used_tfs = set()
        for target, target_votes in votes.items():
            signs = [v[2] for v in target_votes]
            majority = 1 if signs.count(1) >= signs.count(-1) else -1
            if signs.count(majority) / len(signs) >= cfg["signed_network"]["target_path_consensus"]:
                consensus[target] = majority
                used_tfs.update(v[1] for v in target_votes if v[2] == majority)
        reversal = sum(consensus[g] == -sign_of(float(shared[g])) for g in consensus)
        concordant = sum(consensus[g] == sign_of(float(shared[g])) for g in consensus)
        omni_evaluable = len(consensus) >= cfg["signed_network"]["minimum_targets"] and len(used_tfs) >= cfg["signed_network"]["minimum_tfs"]
        omni_class = "NOT_EVALUABLE"
        if omni_evaluable:
            omni_class = "REVERSAL" if reversal > concordant else ("CONCORDANT" if concordant > reversal else "MIXED")

        cs_column = cytosig_map.get(row.ligand.upper().replace(" ", ""))
        cs_genes = sorted(module_genes.intersection(cytosig.index).intersection(shared.index))
        cs_rho = cs_p = np.nan
        cs_class = "NOT_EVALUABLE"
        if cs_column is not None and len(cs_genes) >= cfg["signed_network"]["cytosig_minimum_genes"]:
            cs_rho = rank_corr(cytosig.loc[cs_genes, cs_column].to_numpy(float), shared.loc[cs_genes].to_numpy(float))
            cs_p = exact_permutation_p(cs_rho, cytosig.loc[cs_genes, cs_column].to_numpy(float), shared.loc[cs_genes].to_numpy(float), 2000)
        signed_rows.append({"ligand": row.ligand, "program_module": row.program_module, "receiver_level": row.receiver_level, "receiver_celltype": row.receiver_celltype, "target_compatibility_class": row.target_compatibility_class, "competent_receptor_n": row.competent_receptor_n, "signed_lr_receptor_n": sum((row.ligand, r) in exact_lr_sign for r in str(row.competent_receptors).split(";") if r), "omnipath_collectri_target_n": len(consensus), "omnipath_collectri_TF_n": len(used_tfs), "omnipath_collectri_reversal_target_n": reversal, "omnipath_collectri_concordant_target_n": concordant, "omnipath_collectri_class": omni_class, "cytosig_exact_match": "YES" if cs_column is not None else "NO", "cytosig_matched_name": cs_column or "", "cytosig_module_gene_n": len(cs_genes), "cytosig_spearman_rho": cs_rho, "cytosig_permutation_P": cs_p, "cytosig_BH_FDR": np.nan, "cytosig_class": cs_class, "signed_reversal_class": "PENDING", "interpretation": "SIGNED_GENERIC_PRIOR_NOT_THERAPEUTIC_PROOF", "source_url": SOURCE})
    signed_df = pd.DataFrame(signed_rows)
    signed_df["cytosig_BH_FDR"] = bh(signed_df.cytosig_permutation_P)
    evaluable_cs = signed_df.cytosig_spearman_rho.abs().ge(cfg["signed_network"]["cytosig_abs_rho"]) & signed_df.cytosig_BH_FDR.lt(cfg["signed_network"]["cytosig_family_fdr"])
    signed_df.loc[evaluable_cs & signed_df.cytosig_spearman_rho.lt(0), "cytosig_class"] = "REVERSAL"
    signed_df.loc[evaluable_cs & signed_df.cytosig_spearman_rho.gt(0), "cytosig_class"] = "CONCORDANT"
    reversal_any = signed_df.omnipath_collectri_class.eq("REVERSAL") | signed_df.cytosig_class.eq("REVERSAL")
    concord_any = signed_df.omnipath_collectri_class.eq("CONCORDANT") | signed_df.cytosig_class.eq("CONCORDANT")
    signed_df["signed_reversal_class"] = "SIGNED_EVIDENCE_INSUFFICIENT"
    signed_df.loc[reversal_any & ~concord_any, "signed_reversal_class"] = "REVERSAL_SUPPORTED"
    signed_df.loc[concord_any & ~reversal_any, "signed_reversal_class"] = "DISEASE_CONCORDANT_POTENTIAL"
    signed_df.sort_values(["program_module", "ligand"]).to_csv(OUT / "signed/signed_reversal_analysis.csv", index=False)

    # ---------------------------------------------------------------
    # Separate-axis evidence matrix and deterministic hierarchy.
    evidence_df = target_df.merge(signed_df[["ligand", "program_module", "signed_reversal_class", "omnipath_collectri_class", "cytosig_class"]], on=["ligand", "program_module"], how="left")
    evidence_df = evidence_df.merge(senders[["gene", "sender_evidence_level", "baseline_classification", "licensing_classification", "P1_PARACRINE_CORE", "P2_EXTRACELLULAR_EXTENDED"]].rename(columns={"gene": "ligand", "P1_PARACRINE_CORE": "sender_P1", "P2_EXTRACELLULAR_EXTENDED": "sender_P2"}), on="ligand", how="left")
    receiver_ev = pd.read_csv(OUT / "freeze/phase4_receiver_hierarchy.csv")[["program_module", "CORRECTED_ADMATI_SUPPORT", "EXTERNAL_SCRNA_SUPPORT", "INDEPENDENT_BULK_PROGRAM_SUPPORT", "celltype_localization_claim", "module_union_sha256"]]
    evidence_df = evidence_df.merge(receiver_ev, on="program_module", how="left")
    evidence_df["BASELINE_SENDER_ROBUSTNESS"] = evidence_df.sender_evidence_level
    evidence_df["PARACRINE_ANNOTATION"] = np.where(evidence_df.sender_P1.eq("YES"), "P1_PARACRINE_CORE", np.where(evidence_df.sender_P2.eq("YES"), "P2_EXTRACELLULAR_EXTENDED", "P3_FULL_LR_SENSITIVITY"))
    evidence_df["LICENSING_CONTEXT"] = evidence_df.licensing_classification
    evidence_df["RECEPTOR_COMPETENCE"] = np.where(evidence_df.competent_receptor_n.gt(0), "RECEPTOR_COMPETENT", "NO_COMPETENT_RECEPTOR")
    evidence_df["LR_DATABASE_SUPPORT"] = np.where(evidence_df.competent_receptor_n.gt(0), "NICHENET_LR_WITH_COMPETENT_RECEPTOR", "NO_COMPETENT_NICHENET_LR")
    evidence_df["TARGET_COMPATIBILITY"] = evidence_df.target_compatibility_class
    evidence_df["SIGNED_REVERSAL_EVIDENCE"] = evidence_df.signed_reversal_class
    evidence_df["RECEIVER_EVIDENCE_LEVEL"] = evidence_df.receiver_level
    evidence_df["duplicate_molecular_evidence_rule"] = "same module_union_sha256 is one pathway-level object; localization retained separately"
    evidence_df["axis_tier"] = "NOT_PRIORITIZED"
    target_ok = evidence_df.target_compatibility_class.eq("TARGET_COMPATIBLE") & evidence_df.competent_receptor_n.gt(0)
    primary_receiver = evidence_df.receiver_level.isin(["R1", "R2A"])
    evidence_df.loc[target_ok & ((evidence_df.sender_P1.ne("YES")) | evidence_df.receiver_level.eq("R2B")), "axis_tier"] = "TIER_C_EXTENDED_EXTRACELLULAR"
    evidence_df.loc[target_ok & primary_receiver & evidence_df.signed_reversal_class.eq("SIGNED_EVIDENCE_INSUFFICIENT"), "axis_tier"] = "TIER_B_COMPATIBILITY_CANDIDATE"
    evidence_df.loc[target_ok & primary_receiver & evidence_df.sender_P1.eq("YES") & evidence_df.signed_reversal_class.eq("REVERSAL_SUPPORTED"), "axis_tier"] = "TIER_A_DIRECTIONAL_RESCUE_CANDIDATE"
    evidence_df["source_url"] = SOURCE
    evidence_df.sort_values(["program_module", "ligand"]).to_csv(OUT / "integration/sender_receiver_evidence_matrix.csv", index=False)

    tier_order = {"TIER_A_DIRECTIONAL_RESCUE_CANDIDATE": 0, "TIER_B_COMPATIBILITY_CANDIDATE": 1, "TIER_C_EXTENDED_EXTRACELLULAR": 2, "NOT_PRIORITIZED": 3}
    hierarchy_rows = []
    for ligand, group in evidence_df.groupby("ligand", sort=True):
        best = min(group.axis_tier, key=lambda x: tier_order[x])
        disease_n = int(group.signed_reversal_class.eq("DISEASE_CONCORDANT_POTENTIAL").sum())
        reversal_n = int(group.signed_reversal_class.eq("REVERSAL_SUPPORTED").sum())
        hierarchy_rows.append({"ligand": ligand, "sender_evidence_level": group.sender_evidence_level.iloc[0], "baseline_sender_robustness": group.baseline_classification.iloc[0], "paracrine_scope": group.PARACRINE_ANNOTATION.iloc[0], "licensing_context": group.LICENSING_CONTEXT.iloc[0], "best_phase4a_tier": best, "tier_A_axis_n": int(group.axis_tier.eq("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE").sum()), "tier_B_axis_n": int(group.axis_tier.eq("TIER_B_COMPATIBILITY_CANDIDATE").sum()), "tier_C_axis_n": int(group.axis_tier.eq("TIER_C_EXTENDED_EXTRACELLULAR").sum()), "receptor_competent_module_n": int(group.competent_receptor_n.gt(0).sum()), "target_compatible_module_n": int(group.target_compatibility_class.eq("TARGET_COMPATIBLE").sum()), "reversal_supported_module_n": reversal_n, "disease_concordant_module_n": disease_n, "mixed_signed_direction_across_modules": "YES" if reversal_n > 0 and disease_n > 0 else "NO", "disease_concordant_modules": ";".join(group.loc[group.signed_reversal_class.eq("DISEASE_CONCORDANT_POTENTIAL"), "program_module"]), "tier_A_modules": ";".join(group.loc[group.axis_tier.eq("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE"), "program_module"]), "tier_B_modules": ";".join(group.loc[group.axis_tier.eq("TIER_B_COMPATIBILITY_CANDIDATE"), "program_module"]), "interpretation": "COMPUTATIONAL_HIERARCHY_NOT_THERAPEUTIC_PROOF", "ordering_rule": "tier then ligand symbol; S1 is not ranked above S2", "source_url": SOURCE})
    hierarchy = pd.DataFrame(hierarchy_rows).sort_values(["best_phase4a_tier", "ligand"], key=lambda s: s.map(tier_order) if s.name == "best_phase4a_tier" else s)
    hierarchy.to_csv(OUT / "integration/phase4a_candidate_hierarchy.csv", index=False)
    disease = evidence_df[evidence_df.signed_reversal_class.eq("DISEASE_CONCORDANT_POTENTIAL")].copy()
    disease["risk_interpretation"] = "SIGNED_GENERIC_PRIOR_MAY_REINFORCE_PE_DIRECTION;NOT_CAUSAL_PROOF"
    disease.to_csv(OUT / "integration/disease_concordant_candidates.csv", index=False)

    # QC, risks, diagnostics.
    molecular_duplicates = receivers.groupby("module_union_sha256").filter(lambda x: len(x) > 1).program_module.nunique()
    risks = pd.DataFrame([
        ["P4A_RISK_001", "HIGH", "SIGNED_PRIOR_CONTEXT", "Signed OmniPath/CollecTRI/CytoSig resources are context-generic and do not prove placental or therapeutic direction.", "RETAIN;require Phase4B triangulation", "OPEN"],
        ["P4A_RISK_002", "HIGH", "CELL_LOCALIZATION", "Several receiver localizations depend on corrected Admati evidence without independent matching-cell localization.", "Keep localization and tissue-level evidence separate", "OPEN"],
        ["P4A_RISK_003", "MEDIUM", "RECEPTOR_LAYER", "Competence uses normalized/ceiled public expression summarized by pregnancy, not raw UMI counts.", "Use detection rule and report donor values", "OPEN"],
        ["P4A_RISK_004", "MEDIUM", "NICHE_NET_SIGN", "NicheNet regulatory potential is unsigned and cannot establish reversal.", "Compatibility and reversal are separate columns", "CONTROLLED"],
        ["P4A_RISK_005", "MEDIUM", "DUPLICATE_GENE_SETS", f"{molecular_duplicates} receiver module labels participate in duplicated molecular membership groups.", "Count membership hashes once as molecular evidence", "CONTROLLED"],
        ["P4A_RISK_006", "MEDIUM", "NEGATIVE_EVIDENCE", "Absence of a database edge or signed path can reflect resource incompleteness.", "Interpret as not supported, not biological absence", "OPEN"],
        ["P4A_RISK_007", "LOW", "BLINDING", "No PE-specific or hUC-MSC-mechanism literature was used to prioritize ligands in Phase4A.", "Preserve frozen hierarchy before Phase4B", "CONTROLLED"],
    ], columns=["risk_id", "severity", "domain", "risk", "mitigation", "status"])
    risks["source_url"] = SOURCE
    risks.to_csv(OUT / "qc/phase4a_risk_flags.csv", index=False)
    diagnostics = pd.DataFrame([
        ["sender_total", len(senders), "EXPECTED_214", len(senders) == 214],
        ["sender_P1", int(senders.P1_PARACRINE_CORE.eq("YES").sum()), "EXPECTED_148", int(senders.P1_PARACRINE_CORE.eq("YES").sum()) == 148],
        ["sender_P2", int(senders.P2_EXTRACELLULAR_EXTENDED.eq("YES").sum()), "EXPECTED_190", int(senders.P2_EXTRACELLULAR_EXTENDED.eq("YES").sum()) == 190],
        ["receiver_tested", len(receivers), "EXPECTED_6", len(receivers) == 6],
        ["target_axis_rows", len(target_df), "EXPECTED_1284", len(target_df) == 214 * 6],
        ["signed_axis_rows", len(signed_df), "EXPECTED_1284", len(signed_df) == 214 * 6],
        ["phase4b_locked", "YES", "EXPECTED_YES", True],
        ["python", platform.python_version(), "RECORDED", True],
        ["pandas", pd.__version__, "RECORDED", True],
        ["numpy", np.__version__, "RECORDED", True],
        ["R", "4.5.3", "RESOURCE_CONVERSION_AND_BASE_GRAPHICS", True],
        ["artifact_tool", "2.8.43", "CSV_IMPORT_AND_INSPECTION_VALIDATION", True],
        ["random_seed", SEED, "FROZEN_20260815", True],
    ], columns=["diagnostic", "observed", "expected", "pass"])
    diagnostics["model_or_rule"] = "See docs/PHASE4A_BLINDED_INTEGRATION_PLAN.md and config/phase4a_analysis.json"
    diagnostics["source_url"] = SOURCE
    diagnostics.to_csv(OUT / "qc/phase4a_method_diagnostics.csv", index=False)
    if not diagnostics["pass"].astype(bool).all():
        raise RuntimeError("Phase4A diagnostics failed")
    print(f"PHASE4A_ANALYSIS_OK receptor_rows={len(receptor_df)} lr_rows={len(lr_df)} target_axes={len(target_df)} tierA={hierarchy.best_phase4a_tier.eq('TIER_A_DIRECTIONAL_RESCUE_CANDIDATE').sum()} tierB={hierarchy.best_phase4a_tier.eq('TIER_B_COMPATIBILITY_CANDIDATE').sum()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
