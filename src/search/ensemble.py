"""
Ensemble utilities for AnDOM 2.0.

Three independent arms:
  1. Sequence  — MMseqs2 PSI vs SCOPe 95%  → sccs + CATH code
  2. Structure — ESMFold + Foldseek CATH50 → CATH code + sccs via crosswalk
  3. Profile   — HHblits + HHsearch PDB70  → PDB hits → sccs + CATH via SIFTS

Voting logic
------------
Each arm casts a vote for the top SCOPe class (letter: a/b/c/d/e/f/g).
  3 votes (all arms agree)  → HIGH confidence  🟢
  2 votes (two arms agree)  → MEDIUM confidence 🟡
  1 vote  (single arm only) → LOW confidence   🟠 / 🔵 / 🟣

Ensemble score = weighted sum of normalised -log10(evalue) across contributing arms.
  w_seq=0.35, w_struct=0.40, w_hh=0.25
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


# ── weights ───────────────────────────────────────────────────────────────────
_W_SEQ    = 0.35
_W_STRUCT = 0.40
_W_HH     = 0.25
_MIN_OVERLAP = 0.5


# ── helpers ───────────────────────────────────────────────────────────────────

def _overlap_fraction(s1: int, e1: int, s2: int, e2: int) -> float:
    try:
        s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)
        overlap = max(0, min(e1, e2) - max(s1, s2))
        shorter = min(e1 - s1, e2 - s2)
        return overlap / shorter if shorter > 0 else 0.0
    except Exception:
        return 0.0


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


def _pdb_from_scope(domain_id: str) -> str:
    try:
        return domain_id[1:5].lower()
    except Exception:
        return ""


def _sccs_class(sccs: str) -> str:
    """Return top-level class letter from sccs string e.g. 'a.1.1.2' -> 'a'"""
    if sccs and sccs != "—" and sccs != "?":
        return sccs.split(".")[0]
    return "?"


# ── public: add confidence scores ────────────────────────────────────────────

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


# ── domain architecture bar: 3 rows ──────────────────────────────────────────

def domain_bar_html(
    df_seq: pd.DataFrame | None,
    df_str: pd.DataFrame | None,
    seq_len: int,
    df_hh: pd.DataFrame | None = None,
) -> str:
    lk = lookup.all_domains()

    rows_html = []

    # Row 1 — sequence arm (SCOPe coloured by class)
    row1 = ""
    if df_seq is not None and len(df_seq) > 0:
        for _, row in df_seq.iterrows():
            cls   = lk.get(row["target"], {}).get("cls", "?")
            color = SCOP_COLORS.get(cls, "#888")
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            sccs  = lk.get(row["target"], {}).get("sccs", "?")
            desc  = lk.get(row["target"], {}).get("desc", "")
            tip   = f"{row['target']} | {sccs} | {desc} | e={float(row['evalue']):.1e}"
            row1 += (
                f'<div title="{tip}" style="position:absolute;top:0;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:100%;'
                f'background:{color};opacity:0.9;border-radius:3px;'
                f'border:1px solid rgba(255,255,255,0.6)"></div>'
            )
    rows_html.append(("🧬 Seq", row1, "#3B82F6"))

    # Row 2 — structure arm (teal shades)
    row2 = ""
    if df_str is not None and len(df_str) > 0:
        for _, row in df_str.iterrows():
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            lddt  = float(row.get("lddt", 0.7))
            alpha = 0.5 + 0.5 * lddt
            tip   = (f"{row['target']} | lDDT={lddt:.2f} | "
                     f"e={float(row['evalue']):.1e} | CATH")
            row2 += (
                f'<div title="{tip}" style="position:absolute;top:0;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:100%;'
                f'background:#0F6E56;opacity:{alpha:.2f};border-radius:3px;'
                f'border:1px solid rgba(255,255,255,0.4)"></div>'
            )
    rows_html.append(("🏗 Struct", row2, "#0F6E56"))

    # Row 3 — HHblits arm (purple)
    row3 = ""
    if df_hh is not None and len(df_hh) > 0:
        for _, row in df_hh.iterrows():
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            prob  = float(row.get("prob", 50)) / 100
            alpha = 0.35 + 0.65 * prob
            tip   = (f"{row['hit_name']} | prob={float(row.get('prob',0)):.1f}% | "
                     f"e={float(row['evalue']):.1e} | PDB70")
            row3 += (
                f'<div title="{tip}" style="position:absolute;top:0;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:100%;'
                f'background:#7F77DD;opacity:{alpha:.2f};border-radius:3px;'
                f'border:1px solid rgba(255,255,255,0.4)"></div>'
            )
    rows_html.append(("🔬 Profile", row3, "#7F77DD"))

    # Build combined HTML
    html = '<div style="font-size:11px;color:#888;margin-bottom:4px">Hover bars for hit details · 1 to {seq_len} aa</div>'.replace("{seq_len}", str(seq_len))
    for label, content, color in rows_html:
        has_content = bool(content)
        bg = "#f0f2f6" if has_content else "#f8f9fa"
        html += (
            f'<div style="display:flex;align-items:center;margin-bottom:4px">'
            f'<span style="width:64px;font-size:11px;color:{color};font-weight:500;flex-shrink:0">{label}</span>'
            f'<div style="position:relative;height:22px;background:{bg};'
            f'border-radius:6px;flex:1;border:1px solid #e0e0e0">'
            f'{content}</div></div>'
        )
    return html


# ── two-arm fusion (legacy, kept for batch) ───────────────────────────────────

def fuse_results(
    df_seq:      pd.DataFrame | None,
    df_str:      pd.DataFrame | None,
    w_seq:       float = _W_SEQ,
    w_struct:    float = _W_STRUCT,
    min_overlap: float = _MIN_OVERLAP,
) -> pd.DataFrame:
    return fuse_results_three(df_seq, df_str, None, w_seq, w_struct, _W_HH, min_overlap)


# ── three-arm voting fusion ───────────────────────────────────────────────────

def fuse_results_three(
    df_seq:      pd.DataFrame | None,
    df_str:      pd.DataFrame | None,
    df_hh:       pd.DataFrame | None,
    w_seq:       float = _W_SEQ,
    w_struct:    float = _W_STRUCT,
    w_hh:        float = _W_HH,
    min_overlap: float = _MIN_OVERLAP,
) -> pd.DataFrame:
    """
    Fuse three arms into one ranked ensemble table with voting confidence.

    Each arm contributes hits. Hits from different arms are matched when
    their query regions overlap ≥50%. Voting:
      - 3 arms agree on SCOPe class → HIGH  (evidence='all_three')
      - 2 arms agree                → MEDIUM (evidence='two_arms')
      - 1 arm only                  → LOW   (evidence='seq_only'|'struct_only'|'hhblits_only')

    Returns DataFrame sorted by ensemble_score descending with columns:
        rank, scope_domain, sccs, cath_domain, cath_code, hh_hit,
        evidence, votes, ensemble_score, seq_evalue, struct_evalue,
        hh_evalue, hh_prob, qstart, qend, lddt
    """
    lk = lookup.all_domains()

    # ── build hit lists ───────────────────────────────────────────────────────
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
                "domain":    str(r["target"]),
                "cath_code": get_cath_code(str(r["target"])),
                "evalue":    float(r["evalue"]),
                "lddt":      float(r.get("lddt", 0)),
                "qstart":    int(r.get("qstart", 0)),
                "qend":      int(r.get("qend",   0)),
            })

    hh_hits: list[dict] = []
    if df_hh is not None and len(df_hh) > 0:
        for _, r in df_hh.iterrows():
            hh_hits.append({
                "domain":    str(r.get("hit_name", r.get("pdb", ""))),
                "pdb":       str(r.get("pdb", "")),
                "cath_code": get_cath_code(str(r.get("pdb", ""))),
                "sccs":      str(r.get("sccs", "—")),
                "evalue":    float(r["evalue"]),
                "prob":      float(r.get("prob", 0)),
                "qstart":    int(r.get("qstart", 0)),
                "qend":      int(r.get("qend",   0)),
            })

    # ── normalise within each arm ─────────────────────────────────────────────
    seq_norm    = _normalise_evalues({h["domain"]: h["evalue"] for h in seq_hits})
    struct_norm = _normalise_evalues({h["domain"]: h["evalue"] for h in struct_hits})
    hh_norm     = _normalise_evalues({h["domain"]: h["evalue"] for h in hh_hits})

    rows: list[dict] = []
    matched_struct: set[str] = set()
    matched_hh:     set[str] = set()

    for sh in seq_hits:
        # find overlapping struct hit
        best_st = None
        for th in struct_hits:
            if _overlap_fraction(sh["qstart"], sh["qend"],
                                 th["qstart"], th["qend"]) >= min_overlap:
                if best_st is None or th["evalue"] < best_st["evalue"]:
                    best_st = th

        # find overlapping HHblits hit
        best_hh = None
        for hh in hh_hits:
            if _overlap_fraction(sh["qstart"], sh["qend"],
                                 hh["qstart"], hh["qend"]) >= min_overlap:
                if best_hh is None or hh["evalue"] < best_hh["evalue"]:
                    best_hh = hh

        # ── voting ────────────────────────────────────────────────────────────
        from db.lookup import cath_code_to_scop_class
        votes = 1
        score = w_seq * seq_norm.get(sh["domain"], 0)

        seq_sccs = sh["sccs"] if sh["sccs"] not in ("—","?","") else "?"
        # Strip ~ prefix (CATH-inferred) for voting comparison
        seq_sccs_clean = seq_sccs.lstrip("~") if seq_sccs != "?" else "?"
        seq_cls  = seq_sccs_clean[0] if seq_sccs_clean != "?" else "?"
        seq_cath = sh["cath_code"] or ""
        seq_pdb  = sh["domain"][1:5].lower() if len(sh["domain"]) >= 5 else ""
        arm_sccs  = [seq_sccs]   # full sccs per arm e.g. "a.1.1.2"
        arm_classes = [seq_cls]  # class letter per arm e.g. "a"
        arm_caths   = [seq_cath]
        arm_pdbs    = [seq_pdb]

        if best_st:
            matched_struct.add(best_st["domain"])
            votes += 1
            score += w_struct * struct_norm.get(best_st["domain"], 0)
            st_cls = cath_code_to_scop_class(best_st["cath_code"]) or "?"
            # Structure arm gives us class letter only (no full sccs from CATH crosswalk)
            arm_sccs.append(st_cls)   # only class-level from struct
            arm_classes.append(st_cls)
            arm_caths.append(best_st["cath_code"] or "")
            arm_pdbs.append(best_st["domain"][:4].lower())

        if best_hh:
            matched_hh.add(best_hh["domain"])
            votes += 1
            score += w_hh * hh_norm.get(best_hh["domain"], 0)
            hh_sccs = best_hh["sccs"] if best_hh["sccs"] not in ("—","?","") else "?"
            hh_sccs_clean = hh_sccs.lstrip("~") if hh_sccs != "?" else "?"
            hh_cls  = hh_sccs_clean[0] if hh_sccs_clean != "?" else "?"
            arm_sccs.append(hh_sccs_clean)   # use clean sccs for consensus
            arm_classes.append(hh_cls)
            arm_caths.append(best_hh.get("cath_code", "") or "")
            arm_pdbs.append(best_hh["pdb"][:4].lower() if best_hh.get("pdb") else "")

        # ── Consensus SCOPe: hierarchical agreement ───────────────────────────
        # Try most specific first: full sccs (e.g. "a.1.1.2")
        # Then superfamily (first 3 parts: "a.1.1"), fold ("a.1"), class ("a")
        # Use the most specific level where ≥2 arms agree
        def _majority(values):
            v = [x for x in values if x and x != "?"]
            if not v: return "?"
            best = max(set(v), key=v.count)
            return best if v.count(best) >= 2 else "?"

        # Full sccs agreement (seq + hh both have full sccs)
        full_sccs_vals = [s for s in arm_sccs if "." in s]
        agreed_sccs = _majority(full_sccs_vals)

        # Superfamily level (a.1.1)
        if agreed_sccs == "?":
            sf_vals = [".".join(s.split(".")[:3]) for s in arm_sccs if s.count(".") >= 2]
            agreed_sccs = _majority(sf_vals)

        # Fold level (a.1)
        if agreed_sccs == "?":
            fold_vals = [".".join(s.split(".")[:2]) for s in arm_sccs if "." in s]
            agreed_sccs = _majority(fold_vals)

        # Class letter fallback
        if agreed_sccs == "?":
            agreed_sccs = _majority(arm_classes)

        # Agreement bonus: +0.2 if ≥2 arms agree at any level
        if agreed_sccs != "?":
            score = min(1.0, score + 0.2)

        most_common_cls = agreed_sccs[0] if agreed_sccs and agreed_sccs != "?" else "?"

        # Consensus CATH
        valid_cath = [c for c in arm_caths if c]
        agreed_cath = max(set(valid_cath), key=valid_cath.count) if valid_cath else "—"

        # Consensus PDB
        agreed_pdb = seq_pdb or (arm_pdbs[1] if len(arm_pdbs) > 1 else "") or ""

        if votes == 3:
            evidence = "all_three"
        elif votes == 2:
            evidence = "two_arms"
        else:
            evidence = "seq_only"

        rows.append({
            "scope_domain":  sh["domain"],
            "sccs":          sh["sccs"],
            "cath_domain":   best_st["domain"] if best_st else "—",
            "cath_code":     agreed_cath,
            "hh_hit":        best_hh["domain"] if best_hh else "—",
            "evidence":      evidence,
            "votes":         votes,
            "arms_classes":  " / ".join(arm_sccs),
            "agreed_sccs":   agreed_sccs,
            "agreed_cath":   agreed_cath,
            "agreed_pdb":    agreed_pdb,
            "ensemble_score":round(score, 4),
            "seq_evalue":    sh["evalue"],
            "struct_evalue": best_st["evalue"] if best_st else None,
            "hh_evalue":     best_hh["evalue"] if best_hh else None,
            "hh_prob":       best_hh["prob"] if best_hh else None,
            "qstart":        sh["qstart"],
            "qend":          sh["qend"],
            "lddt":          best_st["lddt"] if best_st else None,
        })

    # struct_only hits
    for th in struct_hits:
        if th["domain"] in matched_struct:
            continue
        # check if any hh overlaps
        best_hh = None
        for hh in hh_hits:
            if _overlap_fraction(th["qstart"], th["qend"],
                                 hh["qstart"], hh["qend"]) >= min_overlap:
                if best_hh is None or hh["evalue"] < best_hh["evalue"]:
                    best_hh = hh
        votes = 2 if best_hh else 1
        if best_hh:
            matched_hh.add(best_hh["domain"])
        score = w_struct * struct_norm.get(th["domain"], 0)
        if best_hh:
            score += w_hh * hh_norm.get(best_hh["domain"], 0)
        from db.lookup import cath_code_to_scop_class
        st_cls = cath_code_to_scop_class(th["cath_code"]) or "?"
        hh_cls = best_hh["sccs"][0] if best_hh and best_hh["sccs"] not in ("—","?","") else "?"
        valid_cls = [c for c in [st_cls, hh_cls] if c != "?"]
        most_common_cls = max(set(valid_cls), key=valid_cls.count) if valid_cls else "?"
        rows.append({
            "scope_domain":  "—",
            "sccs":          "—",
            "cath_domain":   th["domain"],
            "cath_code":     th["cath_code"],
            "hh_hit":        best_hh["domain"] if best_hh else "—",
            "evidence":      "two_arms" if best_hh else "struct_only",
            "votes":         votes,
            "arms_classes":  "/".join([st_cls, hh_cls] if best_hh else [st_cls]),
            "agreed_sccs":   agreed_sccs,
            "agreed_cath":   th["cath_code"] or "—",
            "agreed_pdb":    th["domain"][:4].lower(),
            "ensemble_score":round(score, 4),
            "seq_evalue":    None,
            "struct_evalue": th["evalue"],
            "hh_evalue":     best_hh["evalue"] if best_hh else None,
            "hh_prob":       best_hh["prob"] if best_hh else None,
            "qstart":        th["qstart"],
            "qend":          th["qend"],
            "lddt":          th["lddt"],
        })

    # hhblits_only hits
    for hh in hh_hits:
        if hh["domain"] in matched_hh:
            continue
        hh_cls = hh["sccs"][0] if hh["sccs"] not in ("—","?","") else "?"
        rows.append({
            "scope_domain":  "—",
            "sccs":          hh["sccs"],
            "cath_domain":   "—",
            "cath_code":     hh["cath_code"],
            "hh_hit":        hh["domain"],
            "evidence":      "hhblits_only",
            "votes":         1,
            "arms_classes":  hh_cls,
            "agreed_sccs":   hh_cls,
            "agreed_cath":   hh["cath_code"] or "—",
            "agreed_pdb":    hh["pdb"][:4].lower() if hh.get("pdb") else "—",
            "ensemble_score":round(w_hh * hh_norm.get(hh["domain"], 0), 4),
            "seq_evalue":    None,
            "struct_evalue": None,
            "hh_evalue":     hh["evalue"],
            "hh_prob":       hh["prob"],
            "qstart":        hh["qstart"],
            "qend":          hh["qend"],
            "lddt":          None,
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("ensemble_score", ascending=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df.reset_index(drop=True)


# ── external enrichment ───────────────────────────────────────────────────────

def enrich_with_pfam(uniprot_id: str, timeout: int = 5) -> list[dict]:
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


def fetch_pdb_function(pdb_id: str, timeout: int = 5) -> dict:
    """Fetch functional info for a PDB entry — UniProt ID, gene, function, Pfam."""
    result = {"uniprot_id": None, "gene": "", "organism": "", "function": "", "pfam": []}
    try:
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/1",
            timeout=timeout,
        )
        if r.status_code != 200:
            return result
        uids = (r.json()
                .get("rcsb_polymer_entity_container_identifiers", {})
                .get("uniprot_ids", []))
        if not uids:
            return result
        uid = uids[0]
        result["uniprot_id"] = uid
        uni  = enrich_with_uniprot(uid, timeout=timeout)
        pfam = enrich_with_pfam(uid, timeout=timeout)
        result.update({
            "gene":     uni.get("gene", ""),
            "organism": uni.get("organism", ""),
            "function": uni.get("function", ""),
            "pfam":     pfam,
        })
    except Exception:
        pass
    return result


# ── benchmark helper ──────────────────────────────────────────────────────────

def benchmark_arms(
    seq_results:    list[dict],
    struct_results: list[dict],
    ground_truth:   dict[str, str],
) -> dict:
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
        seq_sccs = [lk.get(h.get("target",""), {}).get("sccs","") for h in sh]

        def _check_seq(sccs_list, arm):
            if sccs_list and sccs_list[0] == true_sccs:
                stats[arm]["rank1"] += 1
            if true_sccs in sccs_list[:5]:
                stats[arm]["top5"] += 1

        _check_seq(seq_sccs, "seq_only")

        from db.lookup import cath_code_to_scop_class
        struct_classes = [cath_code_to_scop_class(get_cath_code(h.get("target",""))) for h in th]
        true_class = true_sccs[0] if true_sccs else ""

        def _check_struct(classes, arm):
            if classes and classes[0] == true_class:
                stats[arm]["rank1"] += 1
            if true_class in classes[:5]:
                stats[arm]["top5"] += 1

        _check_struct(struct_classes, "struct_only")

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
