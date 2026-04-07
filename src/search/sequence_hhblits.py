"""
HHblits/HHsearch-based twilight zone search.
Step 1: query -> HHblits vs UniClust30 -> query.a3m (MSA)
Step 2: query.a3m -> HHsearch vs PDB70 -> remote homologs with SCOPe annotation
Finds homologs at 10-15% sequence identity.
"""
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
    """
    HHblits search: query -> UniClust30 profile -> PDB70 hits.
    Returns DataFrame with hits annotated with SCOPe classification.
    """
    tmp = Path(tmp_dir) if tmp_dir else Path(TMP_DIR) / uuid.uuid4().hex
    tmp.mkdir(parents=True, exist_ok=True)

    query_fa  = tmp / "query.fa"
    query_a3m = tmp / "query.a3m"
    hhr_file  = tmp / "result.hhr"

    if not fasta.strip().startswith(">"):
        fasta = f">query\n{fasta}"
    query_fa.write_text(fasta)

    # Step 1: build profile with HHblits
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

    # Parse .hhr results
    hits = _parse_hhr(hhr_file)
    if not hits:
        return None, None

    df = pd.DataFrame(hits)
    lk = lookup.all_domains()
    # Map PDB code to SCOPe domain
    df["sccs"]       = df["pdb"].map(lambda x: _pdb_to_sccs(x, lk))
    df["class_name"] = df["sccs"].map(lambda x: x.split(".")[0] if x and x != "?" else "?")
    return df, None

def _parse_hhr(hhr_file: Path) -> list:
    """Parse HHsearch .hhr result file."""
    hits = []
    in_hits = False
    with open(hhr_file) as f:
        for line in f:
            if line.startswith(" No Hit"):
                in_hits = True
                continue
            if in_hits:
                if line.strip() == "" or line.startswith(">"):
                    break
                parts = line.split()
                if len(parts) >= 9 and parts[0].isdigit():
                    try:
                        hits.append({
                            "rank":     int(parts[0]),
                            "pdb":      parts[1][:4].lower(),
                            "prob":     float(parts[2]),
                            "evalue":   float(parts[3]),
                            "score":    float(parts[5]),
                            "qstart":   int(parts[7].split("-")[0]),
                            "qend":     int(parts[7].split("-")[1]),
                            "tstart":   int(parts[8].split("-")[0]),
                            "tend":     int(parts[8].split("-")[1]),
                            "PDB link": f"https://www.rcsb.org/structure/{parts[1][:4].lower()}",
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
