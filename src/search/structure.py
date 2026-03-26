"""
Structure-based structural domain assignment.
ESMFold v1 (Meta AI) predicts structure, Foldseek searches against CATH50.

Public API:
    predict_structure(sequence) -> (pdb_path | None, error | None)
    search_cath(pdb_path, evalue) -> (DataFrame | None, error | None)
    run(sequence, evalue) -> (DataFrame | None, error | None)
"""
import os
import subprocess
import requests
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from config import FOLDSEEK, CATH_DB, TMP_DIR, ESMFOLD_API, ESMFOLD_MAXLEN
import db.lookup as lookup


def predict_structure(sequence: str) -> tuple:
    """
    Predict 3D structure via ESMFold public API.

    Parameters
    ----------
    sequence : raw amino acid sequence (no FASTA header, no spaces)

    Returns
    -------
    ("query_struct.pdb", None)  on success
    (None, error_str)           on failure or length exceeded
    """
    clean = "".join(
        l.strip() for l in sequence.splitlines() if not l.startswith(">")
    ).replace(" ", "")

    if len(clean) > ESMFOLD_MAXLEN:
        return None, (
            f"Sequence length {len(clean)} aa exceeds ESMFold public API "
            f"limit ({ESMFOLD_MAXLEN} aa). Disable structural search or "
            f"install ESMFold locally for longer sequences."
        )
    try:
        r = requests.post(
            ESMFOLD_API,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=clean,
            timeout=90,
        )
        if not r.text.startswith("HEADER"):
            return None, f"ESMFold API error: {r.text[:300]}"
        Path("query_struct.pdb").write_text(r.text)
        return "query_struct.pdb", None
    except requests.exceptions.Timeout:
        return None, "ESMFold API timeout (>90s). Try again or use a shorter sequence."
    except Exception as e:
        return None, str(e)


def search_cath(pdb_path: str, evalue: float = 0.001) -> tuple:
    """
    Search a PDB structure against CATH50 using Foldseek.

    Returns
    -------
    (DataFrame, None)  on success — columns: query, target, evalue, bits,
                        lddt, qstart, qend, PDB link
    (None, error_str)  on failure
    (None, None)       if no hits
    """
    out = "struct_results.tsv"
    cmd = (
        f"{FOLDSEEK} easy-search {pdb_path} {CATH_DB} {out} {TMP_DIR} "
        f"--format-output query,target,evalue,bits,lddt,qstart,qend "
        f"-e {evalue} -v 0"
    )
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        return None, f"Foldseek error: {r.stderr[:400]}"
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        return None, None

    df = pd.read_csv(
        out, sep="\t",
        names=["query", "target", "evalue", "bits", "lddt", "qstart", "qend"]
    )
    df["PDB link"] = df["target"].map(lookup.cath_pdb_url)
    return df.sort_values("evalue").reset_index(drop=True), None


def run(sequence: str, evalue: float = 0.001) -> tuple:
    """
    Full structural search pipeline: ESMFold -> Foldseek -> CATH50.
    Returns (DataFrame | None, error | None).
    """
    pdb_path, err = predict_structure(sequence)
    if err:
        return None, err
    return search_cath(pdb_path, evalue)
