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

    # Annotate each hit with SCOPe sccs and CATH code via fast PDB lookups
    df["sccs"]       = df["pdb"].map(lambda x: lookup.pdb_to_sccs(x) or "—")
    df["class_name"] = df["sccs"].map(
        lambda x: SCOP_CLASS_NAMES.get(x[0], "—") if x and x not in ("—", "?") else "—"
    )
    df["cath_code"]  = df["pdb"].map(lambda x: lookup.pdb_to_cath_code(x) or "—")

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
