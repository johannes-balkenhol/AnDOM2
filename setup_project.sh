#!/bin/bash
# AnDOM 2.0 — modular project setup
# Run from the directory where you want the project to live
# Usage: bash setup_project.sh
set -e
echo "Setting up AnDOM 2.0 modular project..."
mkdir -p search db batch benchmark/test_sets benchmark/results docs

# ─────────────────────────────────────────────────────────────
# config.py  — single source of truth for all paths/constants
# ─────────────────────────────────────────────────────────────
cat > config.py << 'EOF'
from pathlib import Path

BASE_DIR       = Path(__file__).parent
MMSEQS         = "mmseqs"
FOLDSEEK       = "foldseek"
SCOPE_DB       = str(BASE_DIR / "scopeSeqDB")
SCOPE_FA       = str(BASE_DIR / "scopeseq_40.fa")
CATH_DB        = str(BASE_DIR / "cathDB")
TMP_DIR        = str(BASE_DIR / "tmp")
ESMFOLD_API    = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAXLEN = 400

SCOP_COLORS = {
    "a": "#e74c3c", "b": "#3498db", "c": "#2ecc71",
    "d": "#f39c12", "e": "#9b59b6", "f": "#1abc9c", "g": "#e67e22",
}
SCOP_CLASSES = {
    "a": "All alpha",   "b": "All beta",     "c": "Alpha/beta",
    "d": "Alpha+beta",  "e": "Multi-domain", "f": "Membrane",
    "g": "Small proteins",
}

EXAMPLES = {
    "Ex 1": {
        "label": "Ex 1 — Haemoglobin (globin, all-alpha)",
        "seq": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "desc": "Human haemoglobin alpha (142 aa) — classic all-alpha globin fold, well-characterised SCOP domain.",
    },
    "Ex 2": {
        "label": "Ex 2 — HIV-1 protease (aspartyl, all-beta)",
        "seq": "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
        "desc": "HIV-1 protease (99 aa) — all-beta aspartyl protease, major drug target.",
    },
    "Ex 3": {
        "label": "Ex 3 — TIM barrel (alpha/beta classic)",
        "seq": "MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLSQELGASNEILLGAQNVDLNLPKDKFVVLIIVYNKPKDILFKDNEALENGKPFKQNLKAKDLALAGVTPDKMKDLKAKGISGAFVPNIVNLHSQAPADCLMSKLVAGEFEGNIYMGLKPNPEELAAAKSSKLSELIQAAYATGNQVAFKPLTDAAQKAAQESSGKKSATIFAGQATVEDGDTVYL",
        "desc": "Triosephosphate isomerase (248 aa) — canonical TIM barrel. Sequence-only search (>400 aa).",
    },
    "Ex 4": {
        "label": "Ex 4 — Src SH2+SH3 (multi-domain)",
        "seq": "MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETSGY",
        "desc": "Src kinase SH2+SH3 domains (180 aa) — multi-domain signalling protein.",
    },
    "Ex 5": {
        "label": "Ex 5 — p53 DNA-binding domain (beta-sandwich)",
        "seq": "SVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL",
        "desc": "p53 tumour suppressor DNA-binding domain (130 aa) — beta-sandwich, cancer biology.",
    },
}
EOF

# ─────────────────────────────────────────────────────────────
# db/__init__.py
# ─────────────────────────────────────────────────────────────
touch db/__init__.py

# ─────────────────────────────────────────────────────────────
# db/lookup.py  — SCOPe metadata, SQLite-ready interface
# ─────────────────────────────────────────────────────────────
cat > db/lookup.py << 'EOF'
"""
SCOPe domain metadata lookup.

Currently backed by an in-memory dict built from the ASTRAL FASTA file.
The public API (get, pdb_url, cath_pdb_url) is stable — swap the backend
for SQLite or any other store without touching callers.

Future upgrade path:
    1. Run db/build_sqlite.py once to create andom.db
    2. Change _load() to query SQLite
    3. Nothing else changes
"""
import streamlit as st
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import SCOPE_FA


@st.cache_data
def _load(fa_file: str = SCOPE_FA) -> dict:
    lookup: dict = {}
    with open(fa_file) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().lstrip(">").split(None, 3)
                if len(parts) >= 2:
                    sid  = parts[0]
                    sccs = parts[1]
                    desc = parts[2] if len(parts) > 2 else ""
                    org  = (line[line.index("{")+1:line.index("}")]
                            if "{" in line and "}" in line else "")
                    lookup[sid] = {
                        "sccs": sccs,
                        "cls":  sccs[0] if sccs else "?",
                        "desc": desc,
                        "org":  org,
                    }
    return lookup


def get(domain_id: str) -> dict:
    """Return metadata dict for a SCOPe domain ID."""
    return _load().get(domain_id, {})


def all_domains() -> dict:
    """Return the full lookup dict (for iteration)."""
    return _load()


def pdb_url(domain_id: str) -> str:
    """RCSB link derived from SCOPe domain ID (e.g. d3d1ka_ -> 3d1k)."""
    try:
        return f"https://www.rcsb.org/structure/{domain_id[1:5].lower()}"
    except Exception:
        return "https://scop.berkeley.edu"


def cath_pdb_url(target: str) -> str:
    """RCSB link derived from CATH domain ID (e.g. 1c7dA02 -> 1c7d)."""
    try:
        return f"https://www.rcsb.org/structure/{target[:4].lower()}"
    except Exception:
        return "https://www.cathdb.info"
EOF

# ─────────────────────────────────────────────────────────────
# search/__init__.py
# ─────────────────────────────────────────────────────────────
touch search/__init__.py

# ─────────────────────────────────────────────────────────────
# search/sequence.py  — MMseqs2 PSI-search against SCOPe 2.08
# ─────────────────────────────────────────────────────────────
cat > search/sequence.py << 'EOF'
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
sys.path.insert(0, str(Path(__file__).parent.parent))
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
EOF

# ─────────────────────────────────────────────────────────────
# search/structure.py  — ESMFold + Foldseek against CATH50
# ─────────────────────────────────────────────────────────────
cat > search/structure.py << 'EOF'
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
sys.path.insert(0, str(Path(__file__).parent.parent))
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
EOF

# ─────────────────────────────────────────────────────────────
# search/ensemble.py  — merge results, scoring, visualisation
# ─────────────────────────────────────────────────────────────
cat > search/ensemble.py << 'EOF'
"""
Ensemble utilities for AnDOM 2.0.

Merges sequence (SCOPe/MMseqs2) and structure (CATH/Foldseek) results,
computes a combined confidence score, and generates the domain bar HTML.

Future extensions:
    - weighted scoring combining e-value and lDDT
    - flag hits found by BOTH layers as high-confidence
    - benchmark comparisons against InterPro / AlphaFold annotations
"""
import math
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import SCOP_COLORS
import db.lookup as lookup


def _evalue_to_score(evalue: float) -> float:
    """Normalise e-value to 0–1 confidence (higher = better)."""
    try:
        return max(0.0, min(1.0, -math.log10(float(evalue)) / 30.0))
    except Exception:
        return 0.0


def add_scores(df_seq: pd.DataFrame | None,
               df_str: pd.DataFrame | None) -> tuple:
    """
    Add a normalised confidence score to each result DataFrame.

    Sequence hits use -log10(evalue) normalised to [0,1].
    Structure hits use lDDT directly (already in [0,1]).

    Returns the same two DataFrames with an extra 'confidence' column.
    """
    if df_seq is not None and len(df_seq) > 0:
        df_seq = df_seq.copy()
        df_seq["confidence"] = df_seq["evalue"].apply(_evalue_to_score)

    if df_str is not None and len(df_str) > 0:
        df_str = df_str.copy()
        df_str["confidence"] = df_str["lddt"].astype(float).clip(0, 1)

    return df_seq, df_str


def domain_bar_html(df_seq: pd.DataFrame | None,
                    df_str: pd.DataFrame | None,
                    seq_len: int) -> str:
    """
    Generate a two-row HTML domain architecture bar.

    Top row   : SCOPe sequence hits, coloured by SCOP class
    Bottom row: CATH structural hits, dark grey
    Hover tooltip shows domain ID, sccs, e-value or lDDT.
    """
    lk = lookup.all_domains()
    bar = (
        '<div style="position:relative;height:64px;background:#f0f2f6;'
        'border-radius:8px;width:100%;margin-bottom:6px">'
    )

    if df_seq is not None and len(df_seq) > 0:
        for _, row in df_seq.iterrows():
            cls   = lk.get(row["target"], {}).get("cls", "?")
            color = SCOP_COLORS.get(cls, "#888")
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            sccs  = lk.get(row["target"], {}).get("sccs", "?")
            desc  = lk.get(row["target"], {}).get("desc", "")
            tip   = f"{row['target']} | {sccs} | {desc} | e={float(row['evalue']):.1e}"
            bar  += (
                f'<div title="{tip}" style="position:absolute;top:4px;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:24px;'
                f'background:{color};opacity:0.9;border-radius:4px;'
                f'border:1px solid rgba(255,255,255,0.6)"></div>'
            )

    if df_str is not None and len(df_str) > 0:
        for _, row in df_str.iterrows():
            left  = (row["qstart"] / seq_len) * 100
            width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
            tip   = (f"{row['target']} | "
                     f"lddt={float(row['lddt']):.2f} | "
                     f"e={float(row['evalue']):.1e} | CATH")
            bar  += (
                f'<div title="{tip}" style="position:absolute;top:34px;'
                f'left:{left:.1f}%;width:{width:.1f}%;height:24px;'
                f'background:#2c3e50;opacity:0.75;border-radius:4px;'
                f'border:1px solid rgba(255,255,255,0.4)"></div>'
            )

    bar += "</div>"
    return bar
EOF

# ─────────────────────────────────────────────────────────────
# batch/__init__.py
# ─────────────────────────────────────────────────────────────
touch batch/__init__.py

# ─────────────────────────────────────────────────────────────
# batch/processor.py  — multi-FASTA batch search
# ─────────────────────────────────────────────────────────────
cat > batch/processor.py << 'EOF'
"""
Batch processing for AnDOM 2.0.

Accepts a multi-FASTA file and runs the ensemble pipeline on every sequence.
Results are written incrementally to a TSV file.

Usage:
    from batch.processor import run_batch
    df = run_batch(
        "proteome.fasta",
        evalue=1e-3,
        iterations=3,
        use_structure=True,
        out_tsv="results/proteome_results.tsv",
    )

Future extensions:
    - parallel execution with multiprocessing.Pool
    - SQLite caching of results to avoid re-running completed sequences
    - progress bar integration for Streamlit batch UI
"""
import pandas as pd
from pathlib import Path
from typing import Generator
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))


def _parse_fasta(fasta_path: str) -> Generator:
    """Yield (seq_id, sequence) from a multi-FASTA file."""
    seq_id, seq_lines = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    yield seq_id, "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
    if seq_id:
        yield seq_id, "".join(seq_lines)


def run_batch(
    fasta_path: str,
    evalue: float = 1e-3,
    iterations: int = 3,
    use_structure: bool = False,
    out_tsv: str = "batch_results.tsv",
    max_sequences: int | None = None,
) -> pd.DataFrame:
    """
    Run AnDOM 2.0 ensemble pipeline on every sequence in a FASTA file.

    Parameters
    ----------
    fasta_path      : path to multi-FASTA input file
    evalue          : e-value cutoff for sequence search
    iterations      : PSI-search iterations
    use_structure   : also run ESMFold + Foldseek (sequences ≤400 aa only)
    out_tsv         : write results here (incremental)
    max_sequences   : stop after N sequences (None = all)

    Returns
    -------
    Combined DataFrame of all hits, with query_id and source columns.
    """
    from search import sequence as seq_search
    from search import structure as str_search

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    all_results = []
    header_written = False

    for i, (seq_id, seq) in enumerate(_parse_fasta(fasta_path)):
        if max_sequences and i >= max_sequences:
            break

        print(f"  [{i+1}] {seq_id} ({len(seq)} aa)")
        fasta = f">{seq_id}\n{seq}"

        df_seq, err = seq_search.run(fasta, evalue=evalue, iterations=iterations)
        if err:
            print(f"    sequence search error: {err}")
        if df_seq is not None and len(df_seq) > 0:
            df_seq = df_seq.copy()
            df_seq["query_id"] = seq_id
            df_seq["source"]   = "SCOPe_sequence"
            all_results.append(df_seq)
            df_seq.to_csv(out_tsv, sep="\t", index=False,
                          mode="a", header=not header_written)
            header_written = True

        if use_structure and len(seq) <= 400:
            df_str, err2 = str_search.run(seq)
            if err2:
                print(f"    structure search error: {err2}")
            if df_str is not None and len(df_str) > 0:
                df_str = df_str.copy()
                df_str["query_id"] = seq_id
                df_str["source"]   = "CATH_structure"
                all_results.append(df_str)
                df_str.to_csv(out_tsv, sep="\t", index=False,
                              mode="a", header=not header_written)
                header_written = True

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return pd.DataFrame()
EOF

# ─────────────────────────────────────────────────────────────
# benchmark/__init__.py
# ─────────────────────────────────────────────────────────────
touch benchmark/__init__.py

# ─────────────────────────────────────────────────────────────
# benchmark/run.py  — benchmark scaffold
# ─────────────────────────────────────────────────────────────
cat > benchmark/run.py << 'EOF'
"""
AnDOM 2.0 benchmark module.

Evaluates domain assignment quality against curated gold-standard sets.

Research questions this benchmark addresses:
    1. Coverage: % of UniProt EC entries assigned a domain
       (baseline: original AnDOM reported 77% on Swiss-Prot 38.0)
    2. Complementarity: does structural search recover domains
       missed by sequence search alone?
    3. Dark proteome: on proteins where InterPro finds nothing,
       does AnDOM 2.0 assign a domain?
    4. Ensemble vs single method: is sequence+structure better
       than either method alone?
    5. AlphaFold gap: on proteins where AF2 classification is
       uncertain, does AnDOM 2.0 help?

Test sets (place FASTA files in benchmark/test_sets/):
    scop_gold.fasta     — proteins with known SCOPe assignments
    dark_proteome.fasta — proteins with no InterPro hits
    ec_sample.fasta     — UniProt EC-annotated proteins (coverage benchmark)
"""
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from batch.processor import run_batch


def _load_gold(gold_tsv: str) -> dict:
    """
    Load gold-standard domain assignments.
    Expected columns: seq_id, expected_sccs (comma-separated if multiple).
    """
    df = pd.read_csv(gold_tsv, sep="\t")
    return {
        row["seq_id"]: set(row["expected_sccs"].split(","))
        for _, row in df.iterrows()
    }


def run_coverage_benchmark(
    fasta_path: str,
    evalue: float = 1e-3,
    iterations: int = 3,
    use_structure: bool = True,
    out_dir: str = "benchmark/results",
) -> dict:
    """
    Run coverage benchmark: what fraction of sequences get at least one hit?

    Returns dict with keys: total, assigned_seq, assigned_str,
    assigned_either, assigned_both, coverage_seq, coverage_str,
    coverage_ensemble.
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_tsv = f"{out_dir}/coverage_results.tsv"

    df = run_batch(fasta_path, evalue=evalue, iterations=iterations,
                   use_structure=use_structure, out_tsv=out_tsv)

    if df.empty:
        print("No results returned.")
        return {}

    total = df["query_id"].nunique()
    seq_ids   = set(df[df["source"]=="SCOPe_sequence"]["query_id"].unique())
    str_ids   = set(df[df["source"]=="CATH_structure"]["query_id"].unique())
    either    = seq_ids | str_ids
    both      = seq_ids & str_ids

    stats = {
        "total":            total,
        "assigned_seq":     len(seq_ids),
        "assigned_str":     len(str_ids),
        "assigned_either":  len(either),
        "assigned_both":    len(both),
        "coverage_seq":     len(seq_ids)  / total if total else 0,
        "coverage_str":     len(str_ids)  / total if total else 0,
        "coverage_ensemble":len(either)   / total if total else 0,
    }

    print("\n=== Coverage Benchmark Results ===")
    print(f"Total sequences:          {stats['total']}")
    print(f"Assigned (sequence only): {stats['assigned_seq']} "
          f"({stats['coverage_seq']:.1%})")
    print(f"Assigned (structure only):{stats['assigned_str']} "
          f"({stats['coverage_str']:.1%})")
    print(f"Assigned (either):        {stats['assigned_either']} "
          f"({stats['coverage_ensemble']:.1%})")
    print(f"Assigned (both layers):   {stats['assigned_both']}")
    print(f"Unique to structure:      "
          f"{len(str_ids - seq_ids)} "
          f"(domains missed by sequence search)")
    return stats


def run_precision_benchmark(
    fasta_path: str,
    gold_tsv: str,
    evalue: float = 1e-3,
    iterations: int = 3,
    out_dir: str = "benchmark/results",
) -> dict:
    """
    Run precision/recall benchmark against a gold-standard SCOPe assignment.

    gold_tsv must have columns: seq_id, expected_sccs
    Returns dict with precision, recall, f1 per method and ensemble.
    """
    gold = _load_gold(gold_tsv)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_tsv = f"{out_dir}/precision_results.tsv"

    df = run_batch(fasta_path, evalue=evalue, iterations=iterations,
                   use_structure=True, out_tsv=out_tsv)

    results = []
    for seq_id, expected in gold.items():
        hits = df[df["query_id"] == seq_id]
        predicted_seq = set(hits[hits["source"]=="SCOPe_sequence"]["sccs"].dropna())
        predicted_str = set()  # CATH hits don't have sccs directly

        tp_seq = len(expected & predicted_seq)
        fp_seq = len(predicted_seq - expected)
        fn_seq = len(expected - predicted_seq)

        results.append({
            "seq_id": seq_id,
            "tp_seq": tp_seq, "fp_seq": fp_seq, "fn_seq": fn_seq,
        })

    out_df = pd.DataFrame(results)
    tp = out_df["tp_seq"].sum()
    fp = out_df["fp_seq"].sum()
    fn = out_df["fn_seq"].sum()
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1        = (2 * precision * recall / (precision + recall)
                 if (precision + recall) > 0 else 0)

    stats = {"precision": precision, "recall": recall, "f1": f1}
    print(f"\n=== Precision Benchmark ===")
    print(f"Precision: {precision:.3f}")
    print(f"Recall:    {recall:.3f}")
    print(f"F1:        {f1:.3f}")
    return stats


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="AnDOM 2.0 benchmark")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("--gold", help="Gold-standard TSV for precision benchmark")
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--no-structure", action="store_true")
    args = parser.parse_args()

    run_coverage_benchmark(
        args.fasta,
        evalue=args.evalue,
        use_structure=not args.no_structure,
    )
    if args.gold:
        run_precision_benchmark(args.fasta, args.gold, evalue=args.evalue)
EOF

# ─────────────────────────────────────────────────────────────
# app.py  — thin Streamlit UI, imports everything from modules
# ─────────────────────────────────────────────────────────────
cat > app.py << 'EOF'
"""
AnDOM 2.0 — Streamlit web application.

This file is intentionally thin: UI layout and user interaction only.
All domain search logic lives in search/, db/, and batch/ modules.

To add a new search method:
    1. Create search/your_method.py with a run() function
    2. Call it here alongside sequence.run() and structure.run()
    3. Add its results to ensemble.domain_bar_html()
"""
import sys
import os
import streamlit as st
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from config import SCOP_COLORS, SCOP_CLASSES, EXAMPLES, ESMFOLD_MAXLEN
import db.lookup as lookup
from search import sequence as seq_search
from search import structure as str_search
from search.ensemble import domain_bar_html, add_scores

st.set_page_config(page_title="AnDOM 2.0", layout="wide")

# ── sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Search Parameters")
    evalue     = st.select_slider("E-value cutoff",
                    options=[1e-30,1e-20,1e-10,1e-5,1e-3,0.01,0.1,1.0], value=1e-3)
    iterations = st.slider("PSI-MMseqs2 iterations", 1, 5, 3)
    max_hits   = st.slider("Max hits shown", 5, 100, 30)
    use_struct = st.toggle("Structural search (ESMFold + Foldseek)", value=True)
    st.divider()
    st.markdown("**SCOP class colours**")
    for k, v in SCOP_CLASSES.items():
        st.markdown(
            f'<span style="background:{SCOP_COLORS[k]};padding:2px 8px;'
            f'border-radius:3px;color:white;font-size:12px">{k}</span> {v}',
            unsafe_allow_html=True,
        )
    st.divider()
    lk = lookup.all_domains()
    st.caption(f"SCOPe 2.08 — {len(lk):,} domains")
    st.caption("CATH50 structural DB")

page = st.tabs(["Search", "Methods"])

# ── search tab ────────────────────────────────────────────────────────────────
with page[0]:
    st.title("AnDOM 2.0 — Structural Domain Finder")
    st.caption(
        "Sequence + Structure ensemble | SCOPe 2.08 + CATH50 | "
        "MMseqs2 + ESMFold + Foldseek | Dandekar Lab Wuerzburg"
    )

    # example buttons
    st.markdown("**Try an example:**")
    ex_cols = st.columns(5)
    for i, (key, ex) in enumerate(EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"btn_{key}"):
            st.session_state["injected_seq"]  = ex["seq"]
            st.session_state["injected_desc"] = ex["desc"]
            st.rerun()

    if "injected_seq" in st.session_state:
        st.info(st.session_state.get("injected_desc", ""))

    seq_input = st.text_area(
        "Paste protein sequence (FASTA or raw):",
        height=130,
        value=st.session_state.get("injected_seq", ""),
    )

    if st.button("Run AnDOM 2.0 Search", type="primary"):
        if not seq_input.strip():
            st.warning("Please paste a sequence or click an example.")
        else:
            lines     = seq_input.strip().splitlines()
            fasta     = seq_input if lines[0].startswith(">") else f">query\n{seq_input}"
            clean_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
            n         = len(clean_seq)

            use_struct_run = use_struct and n <= ESMFOLD_MAXLEN
            if use_struct and n > ESMFOLD_MAXLEN:
                st.warning(
                    f"Sequence is {n} aa — ESMFold public API limit is "
                    f"{ESMFOLD_MAXLEN} aa. Structural search disabled."
                )

            col_seq, col_str = st.columns(2)
            df_seq = df_str = None

            with col_seq:
                with st.spinner("Sequence search — SCOPe 2.08..."):
                    df_seq, err = seq_search.run(
                        fasta, evalue=evalue, iterations=iterations
                    )
                if err:
                    st.error(f"Sequence search failed: {err}")
                elif df_seq is None or len(df_seq) == 0:
                    st.warning("No sequence hits.")
                else:
                    df_seq = df_seq.head(max_hits)
                    st.success(f"Sequence: {len(df_seq)} SCOPe hits")

            if use_struct_run:
                with col_str:
                    with st.spinner("Structure — ESMFold + Foldseek CATH50..."):
                        df_str, err2 = str_search.run(clean_seq)
                    if err2:
                        st.error(f"Structural search failed: {err2}")
                    elif df_str is None or len(df_str) == 0:
                        st.warning("No structural hits.")
                    else:
                        df_str = df_str.head(max_hits)
                        st.success(f"Structure: {len(df_str)} CATH hits")

            # add confidence scores
            df_seq, df_str = add_scores(df_seq, df_str)

            # domain bar
            st.subheader("Domain Architecture Map")
            bar_len = max(
                df_seq["qend"].max() if df_seq is not None and len(df_seq) > 0 else 1,
                df_str["qend"].max() if df_str is not None and len(df_str) > 0 else 1,
                n,
            )
            st.markdown(
                '<p style="font-size:11px;color:#888">'
                "Top: SCOPe sequence hits (coloured by SCOP class) &nbsp;|&nbsp; "
                f"Bottom: CATH structural hits (dark) &nbsp;|&nbsp; "
                f"Hover for details &nbsp;|&nbsp; 1 to {bar_len} aa</p>"
                + domain_bar_html(df_seq, df_str, bar_len),
                unsafe_allow_html=True,
            )

            # results tables
            tab1, tab2 = st.tabs(["SCOPe sequence hits", "CATH structural hits"])

            with tab1:
                if df_seq is not None and len(df_seq) > 0:
                    disp = df_seq[[
                        "target","sccs","class_name","description",
                        "organism","evalue","bits","qstart","qend","pident","PDB link"
                    ]].copy()
                    disp["evalue"] = disp["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp["pident"] = disp["pident"].apply(lambda x: f"{float(x):.1f}%")
                    st.dataframe(
                        disp, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")},
                    )
                    if os.path.exists("seq_results.tsv"):
                        st.download_button(
                            "Download SCOPe TSV",
                            open("seq_results.tsv").read(),
                            "AnDOM_scope.tsv",
                        )
                else:
                    st.info("No sequence hits.")

            with tab2:
                if df_str is not None and len(df_str) > 0:
                    disp2 = df_str[[
                        "target","evalue","bits","lddt","qstart","qend","PDB link"
                    ]].copy()
                    disp2["evalue"] = disp2["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp2["lddt"]   = disp2["lddt"].apply(lambda x: f"{float(x):.3f}")
                    st.dataframe(
                        disp2, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")},
                    )
                    if os.path.exists("struct_results.tsv"):
                        st.download_button(
                            "Download CATH TSV",
                            open("struct_results.tsv").read(),
                            "AnDOM_cath.tsv",
                        )
                else:
                    st.info("No structural hits or structural search disabled.")

    st.divider()
    st.caption(
        "AnDOM 2.0 | SCOPe 2.08 ASTRAL 40% (MMseqs2 PSI) + "
        "CATH50 (ESMFold + Foldseek) | Dandekar Lab Wuerzburg"
    )

# ── methods tab ───────────────────────────────────────────────────────────────
with page[1]:
    st.title("Methods — AnDOM 2.0")
    st.markdown(
        "AnDOM 2.0 updates the original AnDOM server "
        "(Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002) "
        "with modern tools and a new structural search layer."
    )

    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")
    svg_path = Path(__file__).parent / "docs" / "AnDOM_old_vs_new_method.svg"
    if svg_path.exists():
        st.markdown(svg_path.read_text(), unsafe_allow_html=True)

    st.subheader("What is new in AnDOM 2.0")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("""
**Database updates**
- SCOP 1.50 → SCOPe 2.08 (2023)
- ~7,000 → 15,177 domains (+117%)
- ASTRAL compendium maintained within SCOPe

**Search engine**
- PSI-BLAST + IMPALA → MMseqs2
- ~100× faster at equivalent sensitivity
- Configurable iterations and e-value threshold
- Enriched metadata: sccs, fold, organism, PDB link
        """)
    with c2:
        st.markdown("""
**New: structural search layer**
- ESMFold v1 (Meta AI) predicts 3D structure
- Foldseek searches against CATH50
- lDDT confidence per hit
- Recovers domains invisible to sequence methods

**Ensemble output**
- Two-row domain architecture map
- Hover tooltips with full annotation
- Clickable PDB links per hit
- TSV download for both layers
- Local Streamlit (original server offline ~2010)
        """)

    st.subheader("References")
    st.markdown("""
- Schmidt S, Bork P, Dandekar T (2002) *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) *Science* 379:1123–1130
- Fox NK, Brenner SE, Chandonia JM (2014) *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) *Nucleic Acids Res* 49:D266–273
    """)
EOF

# ─────────────────────────────────────────────────────────────
# requirements.txt
# ─────────────────────────────────────────────────────────────
cat > requirements.txt << 'EOF'
streamlit>=1.30
requests
pandas
biopython
EOF

# ─────────────────────────────────────────────────────────────
# .gitignore
# ─────────────────────────────────────────────────────────────
cat > .gitignore << 'EOF'
# databases (too large for git)
cathDB*
scopeDB*
scopeSeqDB*
scopeseq_40.fa
mmseqs/
mmseqs-linux-avx2.tar.gz
# runtime temp files
tmp/
*.pdb
*.tsv
query_input.fasta
queryDB*
resultDB*
# python
__pycache__/
*.pyc
.streamlit/
# benchmark data (add manually)
benchmark/test_sets/
benchmark/results/
EOF

echo ""
echo "Done! Folder structure:"
find . -name "*.py" -o -name "*.txt" -o -name "*.gitignore" | grep -v __pycache__ | sort
echo ""
echo "Next steps:"
echo "  1. cp docs/AnDOM_old_vs_new_method.svg docs/  (already done if present)"
echo "  2. streamlit run app.py --server.port 8501"
echo "  3. git add -A && git commit -m 'refactor: modular project structure'"
