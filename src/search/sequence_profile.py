"""
MMseqs2 two-step profile search against SCOPe 2.08 ASTRAL 95%.
Step 1: Search query against UniRef50 to get sequence profile
Step 2: Search profile against SCOPe 95% DB
Finds homologs at 15-20% identity (twilight zone).
"""
import subprocess
import uuid
from pathlib import Path
import pandas as pd
from config import MMSEQS, SCOPE_DB_95, DATA_DIR, TMP_DIR, SCOP_CLASSES
import db.lookup as lookup

UNIREF50_DB = str(DATA_DIR / "uniref50" / "uniref50")

def run_profile(
    fasta:      str,
    evalue:     float = 0.001,
    iterations: int   = 3,
    threads:    int   = 8,
    tmp_dir:    str | None = None,
) -> tuple:
    """
    Two-step MMseqs2 profile search against SCOPe 95%.
    1. query -> UniRef50 -> profile
    2. profile -> SCOPe 95% -> hits
    """
    tmp = Path(tmp_dir) if tmp_dir else Path(TMP_DIR) / uuid.uuid4().hex
    tmp.mkdir(parents=True, exist_ok=True)

    query_fa   = tmp / "query.fa"
    query_db   = tmp / "queryDB"
    uniref_res = tmp / "uniref_res"
    profile_db = tmp / "profileDB"
    result_db  = tmp / "resultDB"
    result_tsv = tmp / "result.tsv"
    mmseqs_tmp = tmp / "mmseqs_tmp"

    # Ensure FASTA header
    if not fasta.strip().startswith(">"):
        fasta = f">query\n{fasta}"
    query_fa.write_text(fasta)

    cmds = [
        f"{MMSEQS} createdb {query_fa} {query_db} -v 0",
        # Step 1: search against UniRef50
        (f"{MMSEQS} search {query_db} {UNIREF50_DB} {uniref_res} {mmseqs_tmp} "
         f"-s 7 --num-iterations {iterations} -e 0.1 --threads {threads} -v 0"),
        # Step 2: build profile
        f"{MMSEQS} result2profile {query_db} {UNIREF50_DB} {uniref_res} {profile_db} --threads {threads} -v 0",
        # Step 3: search profile against SCOPe 95%
        (f"{MMSEQS} search {profile_db} {SCOPE_DB_95} {result_db} {mmseqs_tmp} "
         f"-e {evalue} --threads {threads} -v 0"),
        # Step 4: convert results
        (f"{MMSEQS} convertalis {profile_db} {SCOPE_DB_95} {result_db} {result_tsv} "
         f"--format-output query,target,evalue,bits,qstart,qend,tstart,tend,pident --threads {threads} -v 0"),
    ]

    for cmd in cmds:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            return None, f"Profile search error: {r.stderr[:400]}"

    if not result_tsv.exists() or result_tsv.stat().st_size == 0:
        return None, None

    df = pd.read_csv(
        result_tsv, sep="\t",
        names=["query","target","evalue","bits","qstart","qend","tstart","tend","pident"]
    )
    lk = lookup.all_domains()
    df["sccs"]        = df["target"].map(lambda x: lk.get(x, {}).get("sccs", "?"))
    df["cls"]         = df["target"].map(lambda x: lk.get(x, {}).get("cls",  "?"))
    df["class_name"]  = df["cls"].map(lambda x: SCOP_CLASSES.get(x, "?"))
    df["description"] = df["target"].map(lambda x: lk.get(x, {}).get("desc", ""))
    df["organism"]    = df["target"].map(lambda x: lk.get(x, {}).get("org",  ""))
    df["PDB link"]    = df["target"].map(lookup.pdb_url)
    return df.sort_values("evalue").reset_index(drop=True), None
