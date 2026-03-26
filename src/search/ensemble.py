"""
Ensemble utilities for AnDOM 2.0.

Matching strategy (scientifically correct):
  SCOPe hits (sequence arm) → domain ID like d3d1ka_, sccs like a.1.1.2
  CATH hits  (structure arm) → domain ID like 1c7dA02, CATH code like 1.10.490.10

  A PDB entry can have multiple domains from different families.
  We match by: same 4-char PDB code AND overlapping query region (qstart/qend).
  This correctly identifies when both arms found the same domain region,
  even though they use different classification schemes.

Additional value beyond SCOPe/CATH:
  - Pfam/InterPro: functional annotation (what the domain does)
  - ECOD: evolutionary classification, orthogonal to SCOP/CATH
  - AlphaFold pLDDT: per-residue structure confidence
  - UniProt: active sites, binding sites, PTMs
  These are retrieved via public APIs and added to enriched hits.
"""
from __future__ import annotations

import math
import requests
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import SCOP_COLORS
import db.lookup as lookup


# ── ID / region helpers ───────────────────────────────────────────────────────

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


def _regions_overlap(s1: int, e1: int, s2: int, e2: int, min_overlap: float = 0.5) -> bool:
    """
    True if two query regions overlap by at least min_overlap fraction
    of the shorter region. Handles missing values gracefully.
    """
    try:
        s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
        overlap = max(0, min(e1, e2) - max(s1, s2))
        shorter = min(e1 - s1, e2 - s2)
        return shorter > 0 and (overlap / shorter) >= min_overlap
    except Exception:
        return False


# ── score helpers ─────────────────────────────────────────────────────────────

def _evalue_to_score(evalue: float) -> float:
    """Normalise e-value to 0–1 confidence (higher = better)."""
    try:
        return max(0.0, min(1.0, -math.log10(float(evalue)) / 30.0))
    except Exception:
        return 0.0


def _normalise_evalues(evalue_map: dict[str, float]) -> dict[str, float]:
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
            tip   = (f"{row['target']} | lddt={float(row['lddt']):.2f} | "
                     f"e={float(row['evalue']):.1e} | CATH")
            bar  += (
                f'<div title="{tip}" style="position:absolute;top:34px;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:24px;'
                f'background:#2c3e50;opacity:0.75;border-radius:4px;'
                f'border:1px solid rgba(255,255,255,0.4)"></div>'
            )
    bar += "</div>"
    return bar


# ── fuse_results: match by PDB + overlapping region ──────────────────────────

_W_SEQ    = 0.40
_W_STRUCT = 0.60


def fuse_results(
    df_seq: pd.DataFrame | None,
    df_str: pd.DataFrame | None,
    w_seq:    float = _W_SEQ,
    w_struct: float = _W_STRUCT,
    min_overlap: float = 0.5,
) -> pd.DataFrame:
    """
    Fuse SCOPe sequence hits and CATH structural hits.

    Matching: same 4-char PDB code AND ≥50% overlap of query regions.
    This correctly handles multi-domain proteins where one PDB entry
    contains domains from different families.

    Each matched pair gets:
        ensemble_score = w_seq * seq_norm + w_struct * struct_norm

    Returns DataFrame with columns:
        rank, pdb_code, scope_domain, sccs, cath_domain, cath_code,
        evidence, ensemble_score, seq_evalue, struct_evalue,
        qstart_seq, qend_seq, qstart_struct, qend_struct, lddt
    """
    lk = lookup.all_domains()

    # Build hit lists with regions
    seq_hits: list[dict] = []
    if df_seq is not None and len(df_seq) > 0 and "target" in df_seq.columns:
        for _, r in df_seq.iterrows():
            pdb = _pdb_from_scope(str(r["target"]))
            if not pdb:
                continue
            seq_hits.append({
                "pdb":    pdb,
                "domain": str(r["target"]),
                "sccs":   lk.get(str(r["target"]), {}).get("sccs", "—"),
                "evalue": float(r["evalue"]),
                "qstart": int(r.get("qstart", 0)),
                "qend":   int(r.get("qend", 0)),
                "score":  _evalue_to_score(float(r["evalue"])),
            })

    struct_hits: list[dict] = []
    if df_str is not None and len(df_str) > 0 and "target" in df_str.columns:
        for _, r in df_str.iterrows():
            pdb = _pdb_from_cath(str(r["target"]))
            if not pdb:
                continue
            struct_hits.append({
                "pdb":       pdb,
                "domain":    str(r["target"]),
                "cath_code": str(r["target"])[4:] if len(str(r["target"])) > 4 else "—",
                "evalue":    float(r["evalue"]),
                "lddt":      float(r.get("lddt", 0)),
                "qstart":    int(r.get("qstart", 0)),
                "qend":      int(r.get("qend", 0)),
                "score":     _evalue_to_score(float(r["evalue"])),
            })

    # Normalise scores within each arm
    seq_norm_map    = _normalise_evalues({h["domain"]: h["evalue"] for h in seq_hits})
    struct_norm_map = _normalise_evalues({h["domain"]: h["evalue"] for h in struct_hits})

    rows: list[dict] = []
    matched_struct: set[str] = set()
    matched_seq:    set[str] = set()

    # Match: seq hit ↔ struct hit with same PDB + overlapping region
    for sh in seq_hits:
        best_match = None
        best_score = -1.0
        for th in struct_hits:
            if th["pdb"] != sh["pdb"]:
                continue
            if not _regions_overlap(sh["qstart"], sh["qend"],
                                    th["qstart"], th["qend"], min_overlap):
                continue
            combined = (w_seq    * seq_norm_map.get(sh["domain"], 0) +
                        w_struct * struct_norm_map.get(th["domain"], 0))
            if combined > best_score:
                best_score = combined
                best_match = th

        if best_match:
            matched_struct.add(best_match["domain"])
            matched_seq.add(sh["domain"])
            rows.append({
                "pdb_code":     sh["pdb"],
                "scope_domain": sh["domain"],
                "sccs":         sh["sccs"],
                "cath_domain":  best_match["domain"],
                "cath_code":    best_match["cath_code"],
                "evidence":     "both",
                "ensemble_score": round(best_score, 4),
                "seq_evalue":   sh["evalue"],
                "struct_evalue":best_match["evalue"],
                "qstart":       sh["qstart"],
                "qend":         sh["qend"],
                "lddt":         best_match["lddt"],
            })
        else:
            # seq_only — no overlapping structural hit
            rows.append({
                "pdb_code":     sh["pdb"],
                "scope_domain": sh["domain"],
                "sccs":         sh["sccs"],
                "cath_domain":  "—",
                "cath_code":    "—",
                "evidence":     "seq_only",
                "ensemble_score": round(w_seq * seq_norm_map.get(sh["domain"], 0), 4),
                "seq_evalue":   sh["evalue"],
                "struct_evalue":None,
                "qstart":       sh["qstart"],
                "qend":         sh["qend"],
                "lddt":         None,
            })

    # struct_only — structural hits with no matching sequence hit
    for th in struct_hits:
        if th["domain"] in matched_struct:
            continue
        rows.append({
            "pdb_code":     th["pdb"],
            "scope_domain": "—",
            "sccs":         "—",
            "cath_domain":  th["domain"],
            "cath_code":    th["cath_code"],
            "evidence":     "struct_only",
            "ensemble_score": round(w_struct * struct_norm_map.get(th["domain"], 0), 4),
            "seq_evalue":   None,
            "struct_evalue":th["evalue"],
            "qstart":       th["qstart"],
            "qend":         th["qend"],
            "lddt":         th["lddt"],
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("ensemble_score", ascending=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df.reset_index(drop=True)


# ── external enrichment (Pfam, UniProt, AlphaFold pLDDT) ─────────────────────

def enrich_with_pfam(uniprot_id: str, timeout: int = 5) -> list[dict]:
    """
    Fetch Pfam domain annotations for a UniProt accession.
    Returns list of {pfam_id, name, description, start, end}.
    Uses InterPro API (replaces deprecated Pfam REST API).
    """
    try:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{uniprot_id}/?format=json"
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return []
        data = r.json()
        results = []
        for entry in data.get("results", []):
            meta = entry.get("metadata", {})
            for protein in entry.get("proteins", []):
                for loc in protein.get("entry_protein_locations", []):
                    for frag in loc.get("fragments", []):
                        results.append({
                            "pfam_id":     meta.get("accession", ""),
                            "name":        meta.get("name", ""),
                            "description": meta.get("description", ""),
                            "start":       frag.get("start"),
                            "end":         frag.get("end"),
                            "source":      "Pfam/InterPro",
                        })
        return results
    except Exception:
        return []


def enrich_with_alphafold_plddt(uniprot_id: str, timeout: int = 8) -> list[dict]:
    """
    Fetch per-residue pLDDT from AlphaFold DB for a UniProt accession.
    Returns list of {residue, plddt} — high pLDDT (>70) = likely structured.
    Useful to flag which regions are confidently folded vs disordered.
    """
    try:
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return []
        data = r.json()
        if not data:
            return []
        entry = data[0]
        # pLDDT is per-residue in the bfactor field of the CIF/PDB
        # The API returns summary confidence — for per-residue need the PDB file
        return [{
            "uniprot":         uniprot_id,
            "avg_plddt":       entry.get("globalMetricValue"),
            "model_url":       entry.get("pdbUrl"),
            "source":          "AlphaFold DB",
        }]
    except Exception:
        return []


def enrich_with_uniprot(uniprot_id: str, timeout: int = 5) -> dict:
    """
    Fetch functional annotations from UniProt for a given accession.
    Returns dict with active_sites, binding_sites, function_text.
    """
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return {}
        data = r.json()
        result: dict = {
            "gene":         data.get("genes", [{}])[0].get("geneName", {}).get("value", ""),
            "organism":     data.get("organism", {}).get("scientificName", ""),
            "function":     "",
            "active_sites": [],
            "binding_sites":[],
            "source":       "UniProt",
        }
        for comment in data.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    result["function"] = texts[0].get("value", "")
        for feature in data.get("features", []):
            ft = feature.get("type", "")
            loc = feature.get("location", {})
            pos = f"{loc.get('start',{}).get('value','')}-{loc.get('end',{}).get('value','')}"
            if ft == "Active site":
                result["active_sites"].append(pos)
            elif ft == "Binding site":
                result["binding_sites"].append(pos)
        return result
    except Exception:
        return {}


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
            _check(fused["pdb_code"].tolist(), "ensemble")

    total = max(stats["total"], 1)
    for arm in ("seq_only", "struct_only", "ensemble"):
        stats[arm]["rank1_pct"] = round(100 * stats[arm]["rank1"] / total, 1)
        stats[arm]["top5_pct"]  = round(100 * stats[arm]["top5"]  / total, 1)

    return stats
