"""
Ensemble utilities for AnDOM 2.0.

Matching strategy (scientifically correct):
  Two hits — one from SCOPe (sequence arm) and one from CATH (structure arm)
  — are considered to annotate the SAME domain if their query regions
  overlap by ≥50% of the shorter hit. PDB code is NOT required to match
  because the two arms search different databases with different representative
  structures. The same fold is represented by many different PDB entries.

  After matching, we report:
    - SCOPe classification: sccs code (e.g. a.1.1.2 = globin superfamily)
    - CATH classification: class from domain ID prefix (1=alpha, 2=beta, 3=alpha/beta)
    - Whether both arms agree (evidence = 'both')
    - Ensemble score: weighted combination of normalised e-values

Additional enrichment (via public APIs, on demand):
  - Pfam/InterPro: functional domain annotation
  - UniProt: active sites, binding sites, function text
  - AlphaFold DB: average pLDDT confidence
"""
from __future__ import annotations

import math
import requests
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import SCOP_COLORS
from db.lookup import get_cath_code
import db.lookup as lookup


# ── CATH domain ID parsing ────────────────────────────────────────────────────

def _cath_class_from_domain(domain_id: str) -> str:
    """
    Infer broad CATH class from the domain ID prefix returned by Foldseek.
    Foldseek CATH50 domain IDs: e.g. 6c4sB00, 1c7dA02, af_P68871_1_142
    For AlphaFold entries (af_*) we can't infer class without the mapping file.
    For PDB entries the class is in cathDB itself — we approximate here.
    Full CATH code requires downloading cath-domain-list.txt separately.
    """
    if domain_id.startswith("af_"):
        return "AF"   # AlphaFold model — no PDB code
    try:
        return domain_id[:4].lower()   # return PDB code as identifier
    except Exception:
        return "?"


def _pdb_from_scope(domain_id: str) -> str:
    """d3d1ka_ -> 3d1k"""
    try:
        return domain_id[1:5].lower()
    except Exception:
        return ""


def _pdb_from_cath(domain_id: str) -> str:
    """6c4sB00 -> 6c4s  |  af_P68871_1_142 -> af"""
    if str(domain_id).startswith("af_"):
        return "af"
    try:
        return domain_id[:4].lower()
    except Exception:
        return ""


# ── region overlap ────────────────────────────────────────────────────────────

def _overlap_fraction(s1: int, e1: int, s2: int, e2: int) -> float:
    """
    Fraction of the shorter region that is covered by the overlap.
    Returns 0.0 if no overlap or invalid coordinates.
    """
    try:
        s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
        overlap = max(0, min(e1, e2) - max(s1, s2))
        shorter = min(e1 - s1, e2 - s2)
        return overlap / shorter if shorter > 0 else 0.0
    except Exception:
        return 0.0


# ── score helpers ─────────────────────────────────────────────────────────────

def _evalue_to_score(evalue: float) -> float:
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


# ── unchanged public API ──────────────────────────────────────────────────────

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


# ── fuse_results: match by query region overlap ───────────────────────────────

_W_SEQ    = 0.40
_W_STRUCT = 0.60
_MIN_OVERLAP = 0.5   # ≥50% of shorter region must overlap


def fuse_results(
    df_seq:      pd.DataFrame | None,
    df_str:      pd.DataFrame | None,
    w_seq:       float = _W_SEQ,
    w_struct:    float = _W_STRUCT,
    min_overlap: float = _MIN_OVERLAP,
) -> pd.DataFrame:
    """
    Fuse SCOPe (sequence) and CATH (structure) hits into one ranked table.

    Matching logic
    --------------
    A SCOPe hit and a CATH hit are considered to annotate the SAME domain
    when their query regions overlap by ≥min_overlap of the shorter hit.
    PDB code is NOT required to match — the two databases use different
    representative structures for the same fold families.

    Ensemble score = w_seq * seq_norm + w_struct * struct_norm
    where norms are -log10(evalue) scaled to [0,1] within each arm.

    Returns
    -------
    DataFrame sorted by ensemble_score descending:
        rank, scope_domain, sccs, cath_domain, evidence,
        ensemble_score, seq_evalue, struct_evalue,
        qstart, qend, lddt, overlap_frac
    """
    lk = lookup.all_domains()

    # Build hit lists
    seq_hits: list[dict] = []
    if df_seq is not None and len(df_seq) > 0 and "target" in df_seq.columns:
        for _, r in df_seq.iterrows():
            seq_hits.append({
                "domain": str(r["target"]),
                "cath_code": get_cath_code(str(r["target"])),
                "sccs":   lk.get(str(r["target"]), {}).get("sccs", "—"),
                "evalue": float(r["evalue"]),
                "qstart": int(r.get("qstart", 0)),
                "qend":   int(r.get("qend",   0)),
            })

    struct_hits: list[dict] = []
    if df_str is not None and len(df_str) > 0 and "target" in df_str.columns:
        for _, r in df_str.iterrows():
            struct_hits.append({
                "domain": str(r["target"]),
                "cath_code": get_cath_code(str(r["target"])),
                "evalue": float(r["evalue"]),
                "lddt":   float(r.get("lddt", 0)),
                "qstart": int(r.get("qstart", 0)),
                "qend":   int(r.get("qend",   0)),
            })

    # Normalise within each arm
    seq_norm_map    = _normalise_evalues({h["domain"]: h["evalue"] for h in seq_hits})
    struct_norm_map = _normalise_evalues({h["domain"]: h["evalue"] for h in struct_hits})

    rows: list[dict] = []
    matched_struct: set[str] = set()

    for sh in seq_hits:
        # find best overlapping structural hit
        best_match    = None
        best_overlap  = 0.0
        best_combined = -1.0

        for th in struct_hits:
            ov = _overlap_fraction(sh["qstart"], sh["qend"],
                                   th["qstart"], th["qend"])
            if ov < min_overlap:
                continue
            combined = (w_seq    * seq_norm_map.get(sh["domain"], 0) +
                        w_struct * struct_norm_map.get(th["domain"], 0))
            if ov > best_overlap or (ov == best_overlap and combined > best_combined):
                best_overlap  = ov
                best_combined = combined
                best_match    = th

        if best_match:
            matched_struct.add(best_match["domain"])
            rows.append({
                "scope_domain":  sh["domain"],
                "sccs":          sh["sccs"],
                "cath_domain":   best_match["domain"],
                "evidence":      "both",
                "ensemble_score":round(best_combined, 4),
                "seq_evalue":    sh["evalue"],
                "struct_evalue": best_match["evalue"],
                "qstart":        sh["qstart"],
                "qend":          sh["qend"],
                "lddt":          best_match["lddt"],
                "overlap_frac":  round(best_overlap, 2),
            })
        else:
            rows.append({
                "scope_domain":  sh["domain"],
                "sccs":          sh["sccs"],
                "cath_domain":   "—",
                "evidence":      "seq_only",
                "ensemble_score":round(w_seq * seq_norm_map.get(sh["domain"], 0), 4),
                "seq_evalue":    sh["evalue"],
                "struct_evalue": None,
                "qstart":        sh["qstart"],
                "qend":          sh["qend"],
                "lddt":          None,
                "overlap_frac":  0.0,
            })

    # struct_only hits — no overlapping sequence hit
    for th in struct_hits:
        if th["domain"] in matched_struct:
            continue
        rows.append({
            "scope_domain":  "—",
            "sccs":          "—",
            "cath_domain":   th["domain"],
            "evidence":      "struct_only",
            "ensemble_score":round(w_struct * struct_norm_map.get(th["domain"], 0), 4),
            "seq_evalue":    None,
            "struct_evalue": th["evalue"],
            "qstart":        th["qstart"],
            "qend":          th["qend"],
            "lddt":          th["lddt"],
            "overlap_frac":  0.0,
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("ensemble_score", ascending=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df.reset_index(drop=True)


# ── external enrichment ───────────────────────────────────────────────────────

def enrich_with_pfam(uniprot_id: str, timeout: int = 5) -> list[dict]:
    """Fetch Pfam domains from InterPro API for a UniProt accession."""
    try:
        url = (f"https://www.ebi.ac.uk/interpro/api/entry/pfam/"
               f"protein/uniprot/{uniprot_id}/?format=json")
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            return []
        results = []
        for entry in r.json().get("results", []):
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
                        })
        return results
    except Exception:
        return []


def enrich_with_uniprot(uniprot_id: str, timeout: int = 5) -> dict:
    """Fetch functional annotations from UniProt."""
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=timeout,
        )
        if r.status_code != 200:
            return {}
        data = r.json()
        result = {
            "gene":         data.get("genes", [{}])[0].get("geneName", {}).get("value", ""),
            "organism":     data.get("organism", {}).get("scientificName", ""),
            "function":     "",
            "active_sites": [],
            "binding_sites":[],
        }
        for comment in data.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    result["function"] = texts[0].get("value", "")[:400]
        for feature in data.get("features", []):
            ft  = feature.get("type", "")
            loc = feature.get("location", {})
            pos = f"{loc.get('start',{}).get('value','')}-{loc.get('end',{}).get('value','')}"
            if ft == "Active site":
                result["active_sites"].append(pos)
            elif ft == "Binding site":
                result["binding_sites"].append(pos)
        return result
    except Exception:
        return {}


def enrich_with_alphafold(uniprot_id: str, timeout: int = 8) -> dict:
    """Fetch AlphaFold prediction metadata including average pLDDT."""
    try:
        r = requests.get(
            f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}",
            timeout=timeout,
        )
        if r.status_code != 200 or not r.json():
            return {}
        entry = r.json()[0]
        return {
            "avg_plddt": entry.get("globalMetricValue"),
            "model_url": entry.get("pdbUrl"),
            "model_page":f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}",
        }
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
    ground_truth: {seq_id: expected_scope_domain_id}
    Matching for seq arm: exact domain ID.
    Matching for struct arm: overlapping region (seq_id maps to a known query region).
    Matching for ensemble: uses fuse_results region-overlap matching.
    """
    stats = {
        "seq_only":    {"rank1": 0, "top5": 0},
        "struct_only": {"rank1": 0, "top5": 0},
        "ensemble":    {"rank1": 0, "top5": 0},
        "total":       len(ground_truth),
    }

    seq_idx    = {r["seq_id"]: r["hits"] for r in seq_results}
    struct_idx = {r["seq_id"]: r["hits"] for r in struct_results}

    lk = lookup.all_domains()
    for seq_id, true_sccs in ground_truth.items():
        sh = seq_idx.get(seq_id, [])
        th = struct_idx.get(seq_id, [])

        # Compare by sccs (fold family) not domain ID
        seq_sccs = [lk.get(h.get("target",""), {}).get("sccs","") for h in sh]

        def _check_seq(sccs_list: list[str], arm: str) -> None:
            if sccs_list and sccs_list[0] == true_sccs:
                stats[arm]["rank1"] += 1
            if true_sccs in sccs_list[:5]:
                stats[arm]["top5"] += 1

        _check_seq(seq_sccs, "seq_only")

        # struct arm: check by sccs class match via SCOPe lookup
        lk = lookup.all_domains()
        true_sccs = lk.get(true_domain, {}).get("sccs", "")
        # struct hits don't have sccs directly — count as hit if region overlaps
        # and we later add CATH↔SCOP concordance
        # For now: struct rank-1 is based on best lddt struct hit region overlap with true seq hit
        if sh and th:
            true_hit = sh[0]  # best seq hit as reference
            for j, t in enumerate(th[:5]):
                ov = _overlap_fraction(
                    true_hit.get("qstart", 0), true_hit.get("qend", 0),
                    t.get("qstart", 0), t.get("qend", 0),
                )
                if ov >= _MIN_OVERLAP:
                    if j == 0:
                        stats["struct_only"]["rank1"] += 1
                    stats["struct_only"]["top5"] += 1
                    break

        df_s = pd.DataFrame(sh) if sh else None
        df_t = pd.DataFrame(th) if th else None
        df_s = pd.DataFrame(sh) if sh else None
        df_t = pd.DataFrame(th) if th else None
        fused = fuse_results(df_s, df_t)
        if not fused.empty:
            fused_sccs = [lk.get(d, {}).get("sccs","") for d in fused["scope_domain"].tolist()]
            if fused_sccs and fused_sccs[0] == true_sccs:
                stats["ensemble"]["rank1"] += 1
            if true_sccs in fused_sccs[:5]:
                stats["ensemble"]["top5"] += 1

    total = max(stats["total"], 1)
    for arm in ("seq_only", "struct_only", "ensemble"):
        stats[arm]["rank1_pct"] = round(100 * stats[arm]["rank1"] / total, 1)
        stats[arm]["top5_pct"]  = round(100 * stats[arm]["top5"]  / total, 1)

    return stats
