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

Test sets (place FASTA files in benchmark/test_sets/):
    scop_gold.fasta     — proteins with known SCOPe assignments
    dark_proteome.fasta — proteins with no InterPro hits
    ec_sample.fasta     — UniProt EC-annotated proteins (coverage benchmark)

New: run_arms_benchmark() for rank-1 / top-5 comparison.
"""
from __future__ import annotations

import sys
import json
import time
import argparse
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from batch.processor import run_batch
from search.ensemble import fuse_results, benchmark_arms


# ── gold-standard helpers ─────────────────────────────────────────────────────

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


def _parse_fasta_headers(fasta_path: str) -> dict[str, str]:
    """
    Parse a SCOPE-style FASTA where the header encodes the classification.
    Example: >d1dlwa_ 1.1.1.1 (A:1-72) ...
    Returns {domain_id: class_or_family}.
    """
    mapping: dict[str, str] = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split()
                if len(parts) >= 2:
                    mapping[parts[0]] = parts[1]
    return mapping


# ── existing benchmarks (unchanged API) ──────────────────────────────────────

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

    total   = df["query_id"].nunique()
    seq_ids = set(df[df["source"] == "SCOPe_sequence"]["query_id"].unique())
    str_ids = set(df[df["source"] == "CATH_structure"]["query_id"].unique())
    either  = seq_ids | str_ids
    both    = seq_ids & str_ids

    stats = {
        "total":             total,
        "assigned_seq":      len(seq_ids),
        "assigned_str":      len(str_ids),
        "assigned_either":   len(either),
        "assigned_both":     len(both),
        "coverage_seq":      len(seq_ids) / total if total else 0,
        "coverage_str":      len(str_ids) / total if total else 0,
        "coverage_ensemble": len(either)  / total if total else 0,
    }

    print("\n=== Coverage Benchmark Results ===")
    print(f"Total sequences:           {stats['total']}")
    print(f"Assigned (sequence only):  {stats['assigned_seq']} "
          f"({stats['coverage_seq']:.1%})")
    print(f"Assigned (structure only): {stats['assigned_str']} "
          f"({stats['coverage_str']:.1%})")
    print(f"Assigned (either):         {stats['assigned_either']} "
          f"({stats['coverage_ensemble']:.1%})")
    print(f"Assigned (both layers):    {stats['assigned_both']}")
    print(f"Unique to structure:       {len(str_ids - seq_ids)} "
          "(domains missed by sequence search)")
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
    Returns dict with precision, recall, f1.
    """
    gold = _load_gold(gold_tsv)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_tsv = f"{out_dir}/precision_results.tsv"

    df = run_batch(fasta_path, evalue=evalue, iterations=iterations,
                   use_structure=True, out_tsv=out_tsv)

    rows = []
    for seq_id, expected in gold.items():
        hits = df[df["query_id"] == seq_id]
        predicted_seq = set(hits[hits["source"] == "SCOPe_sequence"]["sccs"].dropna())
        tp_seq = len(expected & predicted_seq)
        fp_seq = len(predicted_seq - expected)
        fn_seq = len(expected - predicted_seq)
        rows.append({"seq_id": seq_id,
                     "tp_seq": tp_seq, "fp_seq": fp_seq, "fn_seq": fn_seq})

    out_df = pd.DataFrame(rows)
    tp = out_df["tp_seq"].sum()
    fp = out_df["fp_seq"].sum()
    fn = out_df["fn_seq"].sum()
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0)

    stats = {"precision": precision, "recall": recall, "f1": f1}
    print(f"\n=== Precision Benchmark ===")
    print(f"Precision: {precision:.3f}")
    print(f"Recall:    {recall:.3f}")
    print(f"F1:        {f1:.3f}")
    return stats


# ── new: per-arm rank benchmark ───────────────────────────────────────────────

def run_arms_benchmark(
    fasta_path: str,
    evalue: float = 1e-3,
    iterations: int = 3,
    limit: int | None = None,
    out_dir: str = "benchmark/results",
    out_json: str | None = None,
) -> dict:
    """
    Compare sequence-only, structure-only, and ensemble arms on SCOPE-40.

    Uses the SCOPE FASTA headers as ground truth (domain_id == expected hit).
    Computes rank-1 accuracy and top-5 sensitivity for each arm.

    Parameters
    ----------
    fasta_path  : path to scopeseq_40.fa (or any SCOPE-style FASTA)
    limit       : evaluate only first N sequences (None = all)
    out_json    : write JSON result here (optional)

    Returns
    -------
    Dict with seq_only, struct_only, ensemble sub-dicts containing
    rank1, top5, rank1_pct, top5_pct, plus total and elapsed_s.
    """
    from search import sequence as seq_search
    from search import structure as str_search

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    ground_truth = _parse_fasta_headers(fasta_path)

    seqs: dict[str, str] = {}
    with open(fasta_path) as f:
        sid, buf = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sid:
                    seqs[sid] = "".join(buf)
                sid = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
        if sid:
            seqs[sid] = "".join(buf)

    if limit:
        items = list(seqs.items())[:limit]
        seqs = dict(items)
        ground_truth = {k: v for k, v in ground_truth.items() if k in seqs}

    print(f"\n=== Arms Benchmark: {len(seqs)} sequences ===")
    t0 = time.time()

    seq_results:    list[dict] = []
    struct_results: list[dict] = []

    for i, (seq_id, seq) in enumerate(seqs.items(), 1):
        print(f"  [{i}/{len(seqs)}] {seq_id}", end="  ")

        fasta = f">{seq_id}\n{seq}"
        sh, err = seq_search.run(fasta, evalue=evalue, iterations=iterations)
        seq_hits = sh.to_dict("records") if sh is not None and len(sh) > 0 else []
        seq_results.append({"seq_id": seq_id, "hits": seq_hits})
        print(f"seq:{len(seq_hits)}", end="  ")

        if len(seq) <= 400:
            th, err2 = str_search.run(seq)
            struct_hits = th.to_dict("records") if th is not None and len(th) > 0 else []
        else:
            struct_hits = []
        struct_results.append({"seq_id": seq_id, "hits": struct_hits})
        print(f"struct:{len(struct_hits)}")

    elapsed = time.time() - t0
    stats = benchmark_arms(seq_results, struct_results, ground_truth)
    stats["elapsed_s"] = round(elapsed, 1)
    stats["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    # print table
    print(f"\n{'Arm':<14} {'Rank-1':>8} {'Top-5':>8}")
    print("-" * 32)
    for arm in ("seq_only", "struct_only", "ensemble"):
        d = stats[arm]
        print(f"{arm:<14} {d['rank1_pct']:>7.1f}%  {d['top5_pct']:>7.1f}%")
    print(f"\nElapsed: {elapsed:.1f}s")

    if out_json:
        Path(out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(out_json).write_text(json.dumps(stats, indent=2))
        print(f"Results written to {out_json}")

    return stats


# ── CLI ────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AnDOM 2.0 benchmark")
    parser.add_argument("fasta",        help="Input FASTA file")
    parser.add_argument("--gold",       help="Gold-standard TSV for precision benchmark")
    parser.add_argument("--arms",       action="store_true",
                        help="Run rank-1/top-5 per-arm benchmark (SCOPE-style FASTA)")
    parser.add_argument("--evalue",     type=float, default=1e-3)
    parser.add_argument("--limit",      type=int,   default=None,
                        help="Evaluate only first N sequences")
    parser.add_argument("--out-json",   default=None,
                        help="Write arms benchmark JSON here")
    parser.add_argument("--no-structure", action="store_true")
    args = parser.parse_args()

    if args.arms:
        run_arms_benchmark(
            args.fasta,
            evalue=args.evalue,
            limit=args.limit,
            out_json=args.out_json,
        )
    else:
        run_coverage_benchmark(
            args.fasta,
            evalue=args.evalue,
            use_structure=not args.no_structure,
        )
        if args.gold:
            run_precision_benchmark(args.fasta, args.gold, evalue=args.evalue)
