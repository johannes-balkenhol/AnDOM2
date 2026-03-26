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
