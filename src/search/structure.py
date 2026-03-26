"""
Structure-based structural domain assignment.
ESMFold v1 (Meta AI) predicts structure, Foldseek searches against CATH50.

Public API:
    run(sequence, evalue, tmp_dir) -> (DataFrame | None, error | None)

All intermediate files go to tmp_dir so the project root stays clean.
"""
import os
import subprocess
import requests
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import FOLDSEEK, CATH_DB, TMP_DIR, ESMFOLD_API, ESMFOLD_MAXLEN
import db.lookup as lookup


def predict_structure(sequence: str, tmp_dir: Path | None = None) -> tuple:
    """
    Predict 3D structure via ESMFold public API.
    PDB file written to tmp_dir/query_struct.pdb (not project root).
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

        out_dir  = tmp_dir or Path(TMP_DIR)
        out_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = out_dir / "query_struct.pdb"
        pdb_path.write_text(r.text)

        # also keep project-root copy for backwards-compatible download
        try:
            Path("query_struct.pdb").write_text(r.text)
        except Exception:
            pass

        return str(pdb_path), None

    except requests.exceptions.Timeout:
        return None, "ESMFold API timeout (>90s). Try again or use a shorter sequence."
    except Exception as e:
        return None, str(e)


def search_cath(
    pdb_path: str,
    evalue:   float = 0.001,
    tmp_dir:  Path | None = None,
) -> tuple:
    """
    Search a PDB structure against CATH50 using Foldseek.
    All output written to tmp_dir.

    Returns
    -------
    (DataFrame, None)  on success — columns:
        query, target, evalue, bits, lddt, qstart, qend, PDB link
    (None, error_str)  on failure
    (None, None)       if no hits
    """
    out_dir  = tmp_dir or Path(TMP_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv  = out_dir / "struct_results.tsv"
    fs_tmp   = out_dir / "foldseek_tmp"
    fs_tmp.mkdir(exist_ok=True)

    cmd = (
        f"{FOLDSEEK} easy-search {pdb_path} {CATH_DB} {out_tsv} {fs_tmp} "
        f"--format-output query,target,evalue,bits,lddt,qstart,qend "
        f"-e {evalue} -v 0"
    )
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        return None, f"Foldseek error: {r.stderr[:400]}"
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return None, None

    df = pd.read_csv(
        out_tsv, sep="\t",
        names=["query", "target", "evalue", "bits", "lddt", "qstart", "qend"]
    )
    df["PDB link"] = df["target"].map(lookup.cath_pdb_url)

    # backwards-compatible copy
    try:
        import shutil
        shutil.copy(out_tsv, "struct_results.tsv")
    except Exception:
        pass

    return df.sort_values("evalue").reset_index(drop=True), None


def run(
    sequence: str,
    evalue:   float = 0.001,
    tmp_dir:  str | None = None,
) -> tuple:
    """
    Full structural search pipeline: ESMFold -> Foldseek -> CATH50.
    Returns (DataFrame | None, error | None).
    """
    work_dir = Path(tmp_dir) if tmp_dir else None
    pdb_path, err = predict_structure(sequence, tmp_dir=work_dir)
    if err:
        return None, err
    return search_cath(pdb_path, evalue=evalue, tmp_dir=work_dir)
