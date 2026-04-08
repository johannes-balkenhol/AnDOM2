"""
HHblits/HHsearch-based twilight zone search.
Step 1: query -> HHblits vs UniClust30 -> query.a3m (MSA)
Step 2: query.a3m -> HHsearch vs PDB70 -> remote homologs with SCOPe annotation
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

def run_hhblits(
    fasta:      str,
    iterations: int = 3,
    threads:    int = 8,
    tmp_dir:    str | None = None,
) -> tuple:
    tmp = Path(tmp_dir) if tmp_dir else Path(TMP_DIR) / uuid.uuid4().hex
    tmp.mkdir(parents=True, exist_ok=True)

    query_fa  = tmp / "query.fa"
    query_a3m = tmp / "query.a3m"
    hhr_file  = tmp / "result.hhr"

    if not fasta.strip().startswith(">"):
        fasta = f">query\n{fasta}"
    query_fa.write_text(fasta)

    cmd1 = (f"{HHBLITS} -i {query_fa} -d {UNICLUST30_DB} "
            f"-oa3m {query_a3m} -n {iterations} -cpu {threads} -v 0")
    r1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    if r1.returncode != 0:
        return None, f"HHblits error: {r1.stderr[:400]}"
    if not query_a3m.exists() or query_a3m.stat().st_size == 0:
        return None, "HHblits produced no alignment"

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
    lk = lookup.all_domains()
    df["sccs"]       = df["pdb"].map(lambda x: _pdb_to_sccs(x, lk))
    df["class_name"] = df["sccs"].map(lambda x: x.split(".")[0] if x and x != "?" else "?")
    return df, None


def _parse_hhr(hhr_file: Path) -> list:
    """
    Parse HHsearch .hhr summary table using right-split strategy.

    Each hit line (fixed-width) looks like:
      1 6NIL_L DNA dC->dU-editing e 100.0 3.4E-47 3.7E-52  277.7   0.0  104    1-104     1-104 (138)

    Strategy: strip the trailing (tlen), then rsplit() to extract the
    numeric columns from the right. The description in the middle is ignored.
    """
    hits = []
    in_hits = False

    # Matches trailing template length e.g. " (138)"
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
            if not stripped:
                continue
            # Must start with a rank number
            if not stripped[0].isdigit():
                continue

            try:
                # Strip template length from end
                tlen_m = tlen_re.search(line)
                if not tlen_m:
                    continue
                line_body = line[:tlen_m.start()].strip()

                # Split from right to get fixed numeric columns:
                # ... Cols Q_start-Q_end T_start-T_end
                # rsplit to peel off T range, Q range, Cols
                parts_r = line_body.rsplit(None, 3)
                # parts_r[-1] = T_start-T_end  e.g. "1-104"
                # parts_r[-2] = Q_start-Q_end  e.g. "1-104"
                # parts_r[-3] = Cols            e.g. "104"
                # parts_r[-4] = remainder with description + leading numbers
                t_range = parts_r[-1]   # "1-104"
                q_range = parts_r[-2]   # "1-104"
                cols    = parts_r[-3]   # "104"
                left    = parts_r[-4]   # "  1 6NIL_L ... 100.0 3.4E-47 3.7E-52  277.7   0.0"

                # Now rsplit the left part to get SS, Score, P-value, E-value, Prob
                left_parts = left.rsplit(None, 5)
                # left_parts[-1] = SS        "0.0"
                # left_parts[-2] = Score     "277.7"
                # left_parts[-3] = P-value   "3.7E-52"
                # left_parts[-4] = E-value   "3.4E-47"
                # left_parts[-5] = Prob      "100.0"
                # left_parts[0]  = "  1 6NIL_L description..."
                ss      = left_parts[-1]
                score   = left_parts[-2]
                pvalue  = left_parts[-3]
                evalue  = left_parts[-4]
                prob    = left_parts[-5]
                head    = left_parts[0]   # "1 6NIL_L description..."

                # Parse rank and hit name from head
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


def _pdb_to_sccs(pdb: str, lk: dict) -> str:
    """Find SCOPe sccs for a PDB code by scanning lookup."""
    pdb = pdb.lower()
    for domain_id, info in lk.items():
        if domain_id[1:5].lower() == pdb:
            return info.get("sccs", "?")
    return "?"
