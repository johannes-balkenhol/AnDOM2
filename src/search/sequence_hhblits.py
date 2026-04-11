"""
HHblits/HHsearch-based twilight zone search.
Step 1: query -> HHblits vs UniClust30 -> query.a3m (MSA)
Step 2: query.a3m -> HHsearch vs PDB70 -> remote homologs with SCOPe + CATH annotation
Finds homologs at 10-15% sequence identity.
"""
import re
import subprocess
import uuid
from pathlib import Path
import pandas as pd
from config import DATA_DIR, TMP_DIR
import db.lookup as lookup

UNICLUST30_DB = str(DATA_DIR / "uniclust30" / "uniclust30_2018_08" / "uniclust30_2018_08")
PDB70_DB      = str(DATA_DIR / "pdb70" / "pdb70")
SCOP_HHM_DB   = str(DATA_DIR / "scop_hhdb" / "scop_hhdb")  # HHsearch vs SCOPe ASTRAL 95%
HHBLITS       = "hhblits"
HHSEARCH      = "hhsearch"

SCOP_CLASS_NAMES = {
    "a": "All alpha",
    "b": "All beta",
    "c": "Alpha/beta",
    "d": "Alpha+beta",
    "e": "Multi-domain",
    "f": "Membrane",
    "g": "Small proteins",
}


def run_hhblits(
    fasta:      str,
    iterations: int = 3,
    threads:    int = 8,
    tmp_dir:    str | None = None,
) -> tuple:
    """
    HHblits search: query -> UniClust30 profile -> PDB70 hits.
    Returns (DataFrame, error_string).
    DataFrame columns: rank, hit_name, pdb, prob, evalue, score,
                       qstart, qend, tstart, tend, sccs, class_name,
                       cath_code, PDB link
    """
    tmp = Path(tmp_dir) if tmp_dir else Path(TMP_DIR) / uuid.uuid4().hex
    tmp.mkdir(parents=True, exist_ok=True)

    query_fa  = tmp / "query.fa"
    query_a3m = tmp / "query.a3m"
    hhr_file  = tmp / "result.hhr"

    if not fasta.strip().startswith(">"):
        fasta = f">query\n{fasta}"
    query_fa.write_text(fasta)

    # Step 1: build MSA profile with HHblits vs UniClust30
    cmd1 = (f"{HHBLITS} -i {query_fa} -d {UNICLUST30_DB} "
            f"-oa3m {query_a3m} -n {iterations} -cpu {threads} -v 0")
    r1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    if r1.returncode != 0:
        return None, f"HHblits error: {r1.stderr[:400]}"
    if not query_a3m.exists() or query_a3m.stat().st_size == 0:
        return None, "HHblits produced no alignment"

    # Step 2: search profile against PDB70
    cmd2 = (f"{HHSEARCH} -i {query_a3m} -d {PDB70_DB} "
            f"-o {hhr_file} -cpu {threads} -v 0 -p 20 -z 250 -b 250")
    r2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
    if r2.returncode != 0:
        return None, f"HHsearch error: {r2.stderr[:400]}"
    if not hhr_file.exists() or hhr_file.stat().st_size == 0:
        return None, None

    hits = _parse_hhr(hhr_file)
    if not hits:
        return None, None

    df = pd.DataFrame(hits)

    # ── Step 3: HHsearch query profile vs SCOPe ASTRAL 95% for direct sccs ──
    # This gives direct SCOPe classification (a.1.1.2 etc.) at twilight zone sensitivity
    scop_hhr = tmp / "scop_result.hhr"
    scop_db  = Path(SCOP_HHM_DB)

    sccs_map: dict[str, str] = {}   # pdb_hit_name -> sccs from SCOPe DB
    if (scop_db.parent / (scop_db.name + "_a3m.ffdata")).exists():
        cmd3 = (f"{HHSEARCH} -i {query_a3m} -d {SCOP_HHM_DB} "
                f"-o {scop_hhr} -cpu {threads} -v 0 -p 20 -z 250 -b 250")
        r3 = subprocess.run(cmd3, shell=True, capture_output=True, text=True)
        if r3.returncode == 0 and scop_hhr.exists():
            scop_hits = _parse_hhr(scop_hhr)
            for sh in scop_hits:
                # SCOPe domain ID is in hit_name: e.g. "d1hbra_"
                # Parse sccs from the domain ID via lookup
                dom = sh.get("hit_name", "")
                sccs = lookup.all_domains().get(dom, {}).get("sccs", "")
                if not sccs and len(dom) >= 5 and dom[0] == "d":
                    # Try looking up by PDB code
                    sccs = lookup.pdb_to_sccs(dom[1:5])
                if sccs:
                    # Store by region — use (qstart, qend) tuple for overlap matching
                    sccs_map[(sh['qstart'], sh['qend'])] = sccs

    def _find_sccs_by_overlap(qs, qe, min_overlap=0.5):
        """Find best sccs from SCOPe HHsearch results by region overlap."""
        best = ""
        best_ov = 0.0
        for (s, e), sccs in sccs_map.items():
            overlap = max(0, min(qe, e) - max(qs, s))
            shorter = min(qe - qs, e - s)
            if shorter > 0:
                ov = overlap / shorter
                if ov >= min_overlap and ov > best_ov:
                    best_ov = ov
                    best = sccs
        return best

    # ── Annotate PDB70 hits with best available SCOPe sccs ──────────────────
    # Priority: 1) direct PDB→SCOPe, 2) SCOPe HHsearch by region, 3) CATH crosswalk
    CATH_NAMES = {"1":"mainly alpha","2":"mainly beta","3":"alpha/beta","4":"few secondary"}

    def _best_sccs(row) -> tuple[str, str, str]:
        pdb = str(row.get("pdb", ""))
        qs  = int(row.get("qstart", 0))
        qe  = int(row.get("qend", 0))

        # 1. Direct PDB→SCOPe lookup
        sccs = lookup.pdb_to_sccs(pdb) or ""
        if sccs:
            return sccs, SCOP_CLASS_NAMES.get(sccs[0], "—"), "direct"

        # 2. SCOPe HHsearch result for overlapping region
        sccs = _find_sccs_by_overlap(qs, qe)
        if sccs:
            return sccs, SCOP_CLASS_NAMES.get(sccs[0], "—"), "hhsearch_scop"

        # 3. CATH crosswalk fallback
        cc = lookup.pdb_to_cath_code(pdb) or ""
        if cc:
            from db.lookup import cath_code_to_scop_class
            cls = cath_code_to_scop_class(cc)
            if cls:
                return f"~{cls}", CATH_NAMES.get(cls, "CATH inferred"), "cath_inferred"

        return "?", "not in SCOPe or CATH", "unknown"

    enriched = df.apply(_best_sccs, axis=1, result_type="expand")
    df["sccs"]        = enriched[0]
    df["class_name"]  = enriched[1]
    df["sccs_source"] = enriched[2]
    df["cath_code"]   = df["pdb"].map(lambda x: lookup.pdb_to_cath_code(x) or "—")

    # Sort by SCOPe annotation quality: direct > hhsearch_scop > cath_inferred > unknown
    # Within each quality tier, keep original rank order (by prob/evalue)
    source_rank = {"direct": 0, "hhsearch_scop": 1, "cath_inferred": 2, "unknown": 3}
    df["_src_rank"] = df["sccs_source"].map(lambda x: source_rank.get(x, 3))
    df = df.sort_values(["_src_rank", df.index.to_series().rename("orig")],
                        ascending=[True, True]).drop(columns=["_src_rank"])
    df = df.reset_index(drop=True)

    return df, None


def _parse_hhr(hhr_file: Path) -> list:
    """
    Parse HHsearch .hhr summary table using right-split strategy.

    Each hit line (fixed-width) looks like:
      1 6NIL_L DNA dC->dU-editing e 100.0 3.4E-47 3.7E-52  277.7   0.0  104    1-104     1-104 (138)

    Strategy: strip trailing (tlen), then rsplit() to peel off numeric
    columns from the right. The variable-width description is ignored.
    """
    hits = []
    in_hits = False
    tlen_re = re.compile(r'\((\d+)\)\s*$')

    with open(hhr_file) as f:
        for line in f:
            if line.startswith(" No Hit"):
                in_hits = True
                continue
            if not in_hits:
                continue
            if line.startswith(">"):
                break
            stripped = line.strip()
            if not stripped or not stripped[0].isdigit():
                continue

            try:
                tlen_m = tlen_re.search(line)
                if not tlen_m:
                    continue
                line_body = line[:tlen_m.start()].strip()

                # Peel off from right: T_range  Q_range  Cols
                parts_r = line_body.rsplit(None, 3)
                t_range = parts_r[-1]
                q_range = parts_r[-2]
                left    = parts_r[-4]

                # Peel off from left remainder: SS  Score  P-value  E-value  Prob
                left_parts = left.rsplit(None, 5)
                evalue  = left_parts[-4]
                prob    = left_parts[-5]
                score   = left_parts[-2]
                head    = left_parts[0]

                head_parts = head.split(None, 2)
                rank     = int(head_parts[0])
                hit_name = head_parts[1]
                pdb_code = hit_name[:4].lower()

                qstart, qend = map(int, q_range.split("-"))
                tstart, tend = map(int, t_range.split("-"))

                hits.append({
                    "rank":      rank,
                    "hit_name":  hit_name,
                    "pdb":       pdb_code,
                    "prob":      float(prob),
                    "evalue":    float(evalue),
                    "score":     float(score),
                    "qstart":    qstart,
                    "qend":      qend,
                    "tstart":    tstart,
                    "tend":      tend,
                    "PDB link":  f"https://www.rcsb.org/structure/{pdb_code}",
                })
            except Exception:
                pass
    return hits
