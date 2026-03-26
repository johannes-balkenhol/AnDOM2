"""
Ensemble utilities for AnDOM 2.0.

Merges sequence (SCOPe/MMseqs2) and structure (CATH/Foldseek) results,
computes a combined confidence score, and generates the domain bar HTML.

Future extensions:
    - weighted scoring combining e-value and lDDT
    - flag hits found by BOTH layers as high-confidence
    - benchmark comparisons against InterPro / AlphaFold annotations
"""
import math
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import SCOP_COLORS
import db.lookup as lookup


def _evalue_to_score(evalue: float) -> float:
    """Normalise e-value to 0–1 confidence (higher = better)."""
    try:
        return max(0.0, min(1.0, -math.log10(float(evalue)) / 30.0))
    except Exception:
        return 0.0


def add_scores(df_seq: pd.DataFrame | None,
               df_str: pd.DataFrame | None) -> tuple:
    """
    Add a normalised confidence score to each result DataFrame.

    Sequence hits use -log10(evalue) normalised to [0,1].
    Structure hits use lDDT directly (already in [0,1]).

    Returns the same two DataFrames with an extra 'confidence' column.
    """
    if df_seq is not None and len(df_seq) > 0:
        df_seq = df_seq.copy()
        df_seq["confidence"] = df_seq["evalue"].apply(_evalue_to_score)

    if df_str is not None and len(df_str) > 0:
        df_str = df_str.copy()
        df_str["confidence"] = df_str["lddt"].astype(float).clip(0, 1)

    return df_seq, df_str


def domain_bar_html(df_seq: pd.DataFrame | None,
                    df_str: pd.DataFrame | None,
                    seq_len: int) -> str:
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
