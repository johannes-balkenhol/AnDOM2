"""
Ensemble utilities for AnDOM 2.0.

Merges sequence (SCOPe/MMseqs2) and structure (CATH/Foldseek) results,
computes a combined confidence score, and generates the domain bar HTML.

New in this version:
    - fuse_results()    — rank-fused DataFrame combining both arms
    - benchmark_arms()  — rank-1 / top-5 comparison across arms
    - evidence column   — "seq_only" | "struct_only" | "both"
"""
from __future__ import annotations

import math
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import SCOP_COLORS
import db.lookup as lookup


# ── score helpers ─────────────────────────────────────────────────────────────

def _evalue_to_score(evalue: float) -> float:
    """Normalise e-value to 0–1 confidence (higher = better)."""
    try:
        return max(0.0, min(1.0, -math.log10(float(evalue)) / 30.0))
    except Exception:
        return 0.0


def _normalise_evalues(evalue_map: dict[str, float]) -> dict[str, float]:
    """
    Convert {id: evalue} → {id: score ∈ [0,1]}.
    score = -log10(evalue) / max(-log10(evalue)) across the set.
    """
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

    Returns the same two DataFrames with an extra 'confidence' column.
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
    Hover tooltip shows domain ID, sccs, e-value or lDDT.
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


# ── new: rank fusion ──────────────────────────────────────────────────────────

# Weights: structure arm weighted slightly higher (more sensitive for dark proteome)
_W_SEQ    = 0.40
_W_STRUCT = 0.60


def fuse_results(
    df_seq: pd.DataFrame | None,
    df_str: pd.DataFrame | None,
    w_seq:    float = _W_SEQ,
    w_struct: float = _W_STRUCT,
) -> pd.DataFrame:
    """
    Fuse sequence and structure hits into a single ranked DataFrame.

    Algorithm
    ---------
    1. Normalise e-values within each arm to [0,1] using log scale.
    2. For each unique domain_id, ensemble_score = w_seq * seq_norm + w_struct * struct_norm.
    3. Add 'evidence' column: seq_only | struct_only | both.
    4. Sort descending by ensemble_score.

    Returns
    -------
    DataFrame with columns: domain_id, evidence, ensemble_score,
    seq_evalue, struct_evalue, seq_confidence, struct_confidence.
    Empty DataFrame if both inputs are None or empty.
    """
    seq_map:    dict[str, float] = {}   # domain_id → evalue
    struct_map: dict[str, float] = {}

    if df_seq is not None and len(df_seq) > 0 and "target" in df_seq.columns:
        for _, r in df_seq.iterrows():
            seq_map[r["target"]] = min(
                float(r["evalue"]),
                seq_map.get(r["target"], 1.0),
            )

    if df_str is not None and len(df_str) > 0 and "target" in df_str.columns:
        for _, r in df_str.iterrows():
            struct_map[r["target"]] = min(
                float(r["evalue"]),
                struct_map.get(r["target"], 1.0),
            )

    if not seq_map and not struct_map:
        return pd.DataFrame()

    seq_norm    = _normalise_evalues(seq_map)
    struct_norm = _normalise_evalues(struct_map)
    all_ids     = set(seq_norm) | set(struct_norm)

    rows = []
    for did in all_ids:
        s  = seq_norm.get(did, 0.0)
        t  = struct_norm.get(did, 0.0)
        ev = "both" if (did in seq_norm and did in struct_norm) \
             else ("seq_only" if did in seq_norm else "struct_only")
        rows.append({
            "domain_id":        did,
            "evidence":         ev,
            "ensemble_score":   round(w_seq * s + w_struct * t, 4),
            "seq_evalue":       seq_map.get(did),
            "struct_evalue":    struct_map.get(did),
            "seq_confidence":   round(s, 4),
            "struct_confidence":round(t, 4),
        })

    df = pd.DataFrame(rows).sort_values("ensemble_score", ascending=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df.reset_index(drop=True)


# ── new: benchmark helper ─────────────────────────────────────────────────────

def benchmark_arms(
    seq_results:    list[dict],      # [{"seq_id": str, "hits": [{"target": str, ...}]}]
    struct_results: list[dict],
    ground_truth:   dict[str, str],  # {seq_id: true_domain_or_sccs}
) -> dict:
    """
    Compare rank-1 and top-5 accuracy across three arms.

    Parameters
    ----------
    seq_results    : per-sequence hit lists from the sequence arm
    struct_results : per-sequence hit lists from the structure arm
    ground_truth   : {seq_id: expected top hit id or sccs}

    Returns
    -------
    Dict with keys: seq_only, struct_only, ensemble — each a sub-dict
    with rank1, top5, rank1_pct, top5_pct.
    """
    stats = {
        "seq_only":    {"rank1": 0, "top5": 0},
        "struct_only": {"rank1": 0, "top5": 0},
        "ensemble":    {"rank1": 0, "top5": 0},
        "total":       len(ground_truth),
    }

    seq_idx    = {r["seq_id"]: r["hits"] for r in seq_results}
    struct_idx = {r["seq_id"]: r["hits"] for r in struct_results}

    for seq_id, true_hit in ground_truth.items():
        sh = seq_idx.get(seq_id, [])
        th = struct_idx.get(seq_id, [])

        # convert hit dicts to DataFrames for fuse_results
        df_s = pd.DataFrame(sh) if sh else None
        df_t = pd.DataFrame(th) if th else None

        def _check(hits: list[dict], arm: str) -> None:
            ids = [h.get("target", "") for h in hits]
            if ids and ids[0] == true_hit:
                stats[arm]["rank1"] += 1
            if true_hit in ids[:5]:
                stats[arm]["top5"] += 1

        _check(sh, "seq_only")
        _check(th, "struct_only")

        fused = fuse_results(df_s, df_t)
        if not fused.empty:
            fused_ids = fused["domain_id"].tolist()
            if fused_ids and fused_ids[0] == true_hit:
                stats["ensemble"]["rank1"] += 1
            if true_hit in fused_ids[:5]:
                stats["ensemble"]["top5"] += 1

    total = max(stats["total"], 1)
    for arm in ("seq_only", "struct_only", "ensemble"):
        stats[arm]["rank1_pct"] = round(100 * stats[arm]["rank1"] / total, 1)
        stats[arm]["top5_pct"]  = round(100 * stats[arm]["top5"]  / total, 1)

    return stats
