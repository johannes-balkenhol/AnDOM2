"""
Ensemble utilities for AnDOM 2.0.

Merges sequence (SCOPe/MMseqs2) and structure (CATH/Foldseek) results,
computes a combined confidence score, and generates the domain bar HTML.

Key design: SCOPe IDs (d3d1ka_) and CATH IDs (1c7dA02) use different
naming schemes but share the 4-char PDB code. fuse_results() matches
on PDB code so both arms can produce 'both' evidence hits.
"""
from __future__ import annotations

import math
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import SCOP_COLORS
import db.lookup as lookup


# ── ID helpers ────────────────────────────────────────────────────────────────

def _pdb_from_scope(domain_id: str) -> str:
    """d3d1ka_ -> 3d1k"""
    try:
        return domain_id[1:5].lower()
    except Exception:
        return ""


def _pdb_from_cath(domain_id: str) -> str:
    """1c7dA02 -> 1c7d"""
    try:
        return domain_id[:4].lower()
    except Exception:
        return ""


# ── score helpers ─────────────────────────────────────────────────────────────

def _evalue_to_score(evalue: float) -> float:
    """Normalise e-value to 0–1 confidence (higher = better)."""
    try:
        return max(0.0, min(1.0, -math.log10(float(evalue)) / 30.0))
    except Exception:
        return 0.0


def _normalise_evalues(evalue_map: dict[str, float]) -> dict[str, float]:
    """Convert {id: evalue} → {id: score ∈ [0,1]} using log scale."""
    if not evalue_map:
        return {}
    EPS = 1e-300
    log_scores = {k: -math.log10(max(v, EPS)) for k, v in evalue_map.items()}
    max_log = max(log_scores.values()) or 1.0
    return {k: v / max_log for k, v in log_scores.items()}


# ── public API (unchanged) ────────────────────────────────────────────────────

def add_scores(
    df_seq: pd.DataFrame | None,
    df_str: pd.DataFrame | None,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """
    Add a normalised confidence score to each result DataFrame.
    Sequence hits  : -log10(evalue) normalised to [0,1]
    Structure hits : lDDT directly (already in [0,1])
    """
    if df_seq is not None and len(df_seq) > 0:
        df_seq = df_seq.copy()
        df_seq["confidence"] = df_seq["evalue"].apply(_evalue_to_score)

    if df_str is not None and len(df_str) > 0:
        df_str = df_str.copy()
        df_str["confidence"] = df_str["lddt"].astype(float).clip(0, 1)

    return df_seq, df_str


def domain_bar_html(
    df_seq: pd.DataFrame | None,
    df_str: pd.DataFrame | None,
    seq_len: int,
) -> str:
    """
    Generate a two-row HTML domain architecture bar.
    Top row   : SCOPe sequence hits, coloured by SCOP class
    Bottom row: CATH structural hits, dark grey
    """
    lk = lookup.all_domains()
    bar = (
        '<div style="position:relative;height:64px;background:#f0f2f6;'
        'border-radius:8px;width:100%;margin-bottom:6px">'
    )

    if df_seq is not None and len(df_seq) > 0:
        for _, row in df_seq.iterrows():
            cls   = lk.get(row["target"], {}).get("cls", "?")
            color = SCOP_COLORS.get(cls, "#888")
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            sccs  = lk.get(row["target"], {}).get("sccs", "?")
            desc  = lk.get(row["target"], {}).get("desc", "")
            tip   = f"{row['target']} | {sccs} | {desc} | e={float(row['evalue']):.1e}"
            bar  += (
                f'<div title="{tip}" style="position:absolute;top:4px;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:24px;'
                f'background:{color};opacity:0.9;border-radius:4px;'
                f'border:1px solid rgba(255,255,255,0.6)"></div>'
            )

    if df_str is not None and len(df_str) > 0:
        for _, row in df_str.iterrows():
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            tip   = (f"{row['target']} | "
                     f"lddt={float(row['lddt']):.2f} | "
                     f"e={float(row['evalue']):.1e} | CATH")
            bar  += (
                f'<div title="{tip}" style="position:absolute;top:34px;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:24px;'
                f'background:#2c3e50;opacity:0.75;border-radius:4px;'
                f'border:1px solid rgba(255,255,255,0.4)"></div>'
            )

    bar += "</div>"
    return bar


# ── fuse_results: match by PDB code ──────────────────────────────────────────

_W_SEQ    = 0.40
_W_STRUCT = 0.60


def fuse_results(
    df_seq: pd.DataFrame | None,
    df_str: pd.DataFrame | None,
    w_seq:    float = _W_SEQ,
    w_struct: float = _W_STRUCT,
) -> pd.DataFrame:
    """
    Fuse sequence (SCOPe) and structure (CATH) hits into a single ranked table.

    Matching strategy
    -----------------
    SCOPe domain IDs (d3d1ka_) and CATH domain IDs (1c7dA02) use different
    naming schemes but share a 4-char PDB code (3d1k / 1c7d).
    We group by PDB code so hits from the same structure are recognised as
    'both' evidence regardless of domain ID format.

    Each PDB code gets:
        ensemble_score = w_seq * seq_norm + w_struct * struct_norm
    where norms are computed within each arm using -log10(evalue).

    Returns
    -------
    DataFrame sorted by ensemble_score descending, with columns:
    rank, pdb_code, scope_domain, cath_domain, evidence,
    ensemble_score, seq_evalue, struct_evalue, lddt.
    """
    # Build {pdb_code: best_evalue} maps for each arm
    seq_map:    dict[str, float] = {}   # pdb -> evalue
    struct_map: dict[str, float] = {}
    seq_meta:   dict[str, dict]  = {}   # pdb -> row data
    struct_meta:dict[str, dict]  = {}

    if df_seq is not None and len(df_seq) > 0 and "target" in df_seq.columns:
        for _, r in df_seq.iterrows():
            pdb = _pdb_from_scope(str(r["target"]))
            if not pdb:
                continue
            ev = float(r["evalue"])
            if pdb not in seq_map or ev < seq_map[pdb]:
                seq_map[pdb]  = ev
                seq_meta[pdb] = r.to_dict()

    if df_str is not None and len(df_str) > 0 and "target" in df_str.columns:
        for _, r in df_str.iterrows():
            pdb = _pdb_from_cath(str(r["target"]))
            if not pdb:
                continue
            ev = float(r["evalue"])
            if pdb not in struct_map or ev < struct_map[pdb]:
                struct_map[pdb]  = ev
                struct_meta[pdb] = r.to_dict()

    if not seq_map and not struct_map:
        return pd.DataFrame()

    seq_norm    = _normalise_evalues(seq_map)
    struct_norm = _normalise_evalues(struct_map)
    all_pdbs    = set(seq_norm) | set(struct_norm)

    rows = []
    for pdb in all_pdbs:
        s  = seq_norm.get(pdb, 0.0)
        t  = struct_norm.get(pdb, 0.0)
        in_seq    = pdb in seq_norm
        in_struct = pdb in struct_norm
        evidence  = "both" if (in_seq and in_struct) else \
                    ("seq_only" if in_seq else "struct_only")

        sm = seq_meta.get(pdb, {})
        tm = struct_meta.get(pdb, {})

        rows.append({
            "pdb_code":       pdb,
            "scope_domain":   sm.get("target", "—"),
            "cath_domain":    tm.get("target", "—"),
            "evidence":       evidence,
            "ensemble_score": round(w_seq * s + w_struct * t, 4),
            "seq_evalue":     seq_map.get(pdb),
            "struct_evalue":  struct_map.get(pdb),
            "lddt":           tm.get("lddt"),
            "sccs":           sm.get("sccs", "—"),
            "seq_norm":       round(s, 4),
            "struct_norm":    round(t, 4),
        })

    df = pd.DataFrame(rows).sort_values("ensemble_score", ascending=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df.reset_index(drop=True)


# ── benchmark helper ──────────────────────────────────────────────────────────

def benchmark_arms(
    seq_results:    list[dict],
    struct_results: list[dict],
    ground_truth:   dict[str, str],
) -> dict:
    """
    Compare rank-1 and top-5 accuracy across three arms.
    ground_truth: {seq_id: expected_pdb_code (4 chars, lowercase)}
    """
    stats = {
        "seq_only":    {"rank1": 0, "top5": 0},
        "struct_only": {"rank1": 0, "top5": 0},
        "ensemble":    {"rank1": 0, "top5": 0},
        "total":       len(ground_truth),
    }

    seq_idx    = {r["seq_id"]: r["hits"] for r in seq_results}
    struct_idx = {r["seq_id"]: r["hits"] for r in struct_results}

    for seq_id, true_pdb in ground_truth.items():
        sh = seq_idx.get(seq_id, [])
        th = struct_idx.get(seq_id, [])

        seq_pdbs    = [_pdb_from_scope(h.get("target", "")) for h in sh]
        struct_pdbs = [_pdb_from_cath(h.get("target", ""))  for h in th]

        def _check(pdbs: list[str], arm: str) -> None:
            if pdbs and pdbs[0] == true_pdb:
                stats[arm]["rank1"] += 1
            if true_pdb in pdbs[:5]:
                stats[arm]["top5"] += 1

        _check(seq_pdbs,    "seq_only")
        _check(struct_pdbs, "struct_only")

        df_s = pd.DataFrame(sh) if sh else None
        df_t = pd.DataFrame(th) if th else None
        fused = fuse_results(df_s, df_t)
        if not fused.empty:
            fused_pdbs = fused["pdb_code"].tolist()
            _check(fused_pdbs, "ensemble")

    total = max(stats["total"], 1)
    for arm in ("seq_only", "struct_only", "ensemble"):
        stats[arm]["rank1_pct"] = round(100 * stats[arm]["rank1"] / total, 1)
        stats[arm]["top5_pct"]  = round(100 * stats[arm]["top5"]  / total, 1)

    return stats
