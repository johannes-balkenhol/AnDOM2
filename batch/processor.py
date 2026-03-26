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
