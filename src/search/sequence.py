"""
Sequence-based structural domain assignment.
MMseqs2 PSI-search against SCOPe 2.08 ASTRAL (40% identity cutoff).

Public API:
    run(fasta, evalue, iterations, tmp_dir) -> (DataFrame | None, error | None)

All intermediate files (queryDB*, resultDB*, seq_results.tsv, query_input.fasta)
are written to tmp_dir so the project root stays clean during multi-user operation.
"""
import os
import subprocess
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import MMSEQS, SCOPE_DB, TMP_DIR, SCOP_CLASSES
import db.lookup as lookup


def run(
    fasta:      str,
    evalue:     float = 1e-3,
    iterations: int   = 3,
    threads:    int   = 8,
    tmp_dir:    str | None = None,
) -> tuple:
    """
    Run MMseqs2 PSI-search against SCOPe.

    Parameters
    ----------
    fasta       : FASTA string (single sequence, with or without header)
    evalue      : E-value threshold
    iterations  : PSI-search iterations (1-5)
    threads     : CPU threads for MMseqs2
    tmp_dir     : directory for all intermediate files (default: TMP_DIR from config)
                  Each job should pass its own isolated directory.

    Returns
    -------
    (DataFrame, None)  on success — DataFrame has columns:
        query, target, evalue, bits, qstart, qend, tstart, tend, pident,
        sccs, cls, class_name, description, organism, PDB link
    (None, error_str)  on failure
    (None, None)       if no hits found
    """
    # ── job-specific working directory ────────────────────────────────────────
    work_dir = Path(tmp_dir) if tmp_dir else Path(TMP_DIR)
    work_dir.mkdir(parents=True, exist_ok=True)

    query_fasta  = work_dir / "query_input.fasta"
    query_db     = work_dir / "queryDB"
    result_db    = work_dir / "resultDB"
    result_tsv   = work_dir / "seq_results.tsv"
    mmseqs_tmp   = work_dir / "mmseqs_tmp"
    mmseqs_tmp.mkdir(exist_ok=True)

    query_fasta.write_text(fasta)

    cmds = [
        f"{MMSEQS} createdb {query_fasta} {query_db} -v 0",
        (f"{MMSEQS} search {query_db} {SCOPE_DB} {result_db} {mmseqs_tmp} "
         f"--num-iterations {iterations} -e {evalue} -v 0 --threads {threads}"),
        (f"{MMSEQS} convertalis {query_db} {SCOPE_DB} {result_db} {result_tsv} "
         f"--format-output query,target,evalue,bits,qstart,qend,tstart,tend,pident -v 0"),
    ]
    for cmd in cmds:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            return None, f"MMseqs2 error: {r.stderr[:400]}"

    if not result_tsv.exists() or result_tsv.stat().st_size == 0:
        return None, None

    df = pd.read_csv(
        result_tsv, sep="\t",
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

    # also write to project root for backwards-compatible download buttons
    try:
        import shutil
        shutil.copy(result_tsv, "seq_results.tsv")
    except Exception:
        pass

    return df.sort_values("evalue").reset_index(drop=True), None
