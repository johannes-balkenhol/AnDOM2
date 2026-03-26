"""
Sequence-based structural domain assignment.
MMseqs2 PSI-search against SCOPe 2.08 ASTRAL (40% identity cutoff).

Public API:
    run(fasta, evalue, iterations) -> (DataFrame | None, error | None)
"""
import os
import glob
import subprocess
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import MMSEQS, SCOPE_DB, TMP_DIR, SCOP_CLASSES
import db.lookup as lookup


def _clean_tmp() -> None:
    """Remove stale MMseqs2 query/result files before each run."""
    for pattern in ("queryDB*", "resultDB*"):
        for f in glob.glob(pattern):
            try:
                os.remove(f)
            except OSError:
                pass
    for f in ("seq_results.tsv", "query_input.fasta"):
        try:
            os.remove(f)
        except OSError:
            pass


def run(fasta: str, evalue: float = 1e-3, iterations: int = 3,
        threads: int = 8) -> tuple:
    """
    Run MMseqs2 PSI-search against SCOPe.

    Parameters
    ----------
    fasta       : FASTA string (single sequence, with or without header)
    evalue      : E-value threshold
    iterations  : PSI-search iterations (1-5)
    threads     : CPU threads for MMseqs2

    Returns
    -------
    (DataFrame, None)  on success — DataFrame has columns:
        query, target, evalue, bits, qstart, qend, tstart, tend, pident,
        sccs, cls, class_name, description, organism, PDB link
    (None, error_str)  on failure
    (None, None)       if no hits found
    """
    _clean_tmp()
    Path("query_input.fasta").write_text(fasta)

    cmds = [
        f"{MMSEQS} createdb query_input.fasta queryDB -v 0",
        (f"{MMSEQS} search queryDB {SCOPE_DB} resultDB {TMP_DIR} "
         f"--num-iterations {iterations} -e {evalue} -v 0 --threads {threads}"),
        (f"{MMSEQS} convertalis queryDB {SCOPE_DB} resultDB seq_results.tsv "
         f"--format-output query,target,evalue,bits,qstart,qend,tstart,tend,pident -v 0"),
    ]
    for cmd in cmds:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            return None, f"MMseqs2 error: {r.stderr[:400]}"

    if not os.path.exists("seq_results.tsv") or os.path.getsize("seq_results.tsv") == 0:
        return None, None

    df = pd.read_csv(
        "seq_results.tsv", sep="\t",
        names=["query", "target", "evalue", "bits",
               "qstart", "qend", "tstart", "tend", "pident"]
    )
    lk = lookup.all_domains()
    df["sccs"]        = df["target"].map(lambda x: lk.get(x, {}).get("sccs", "?"))
    df["cls"]         = df["target"].map(lambda x: lk.get(x, {}).get("cls",  "?"))
    df["class_name"]  = df["cls"].map(lambda x: SCOP_CLASSES.get(x, "?"))
    df["description"] = df["target"].map(lambda x: lk.get(x, {}).get("desc", ""))
    df["organism"]    = df["target"].map(lambda x: lk.get(x, {}).get("org",  ""))
    df["PDB link"]    = df["target"].map(lookup.pdb_url)

    return df.sort_values("evalue").reset_index(drop=True), None
