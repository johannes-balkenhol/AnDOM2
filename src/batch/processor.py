"""
Batch processing for AnDOM 2.0.

Accepts a multi-FASTA file or raw FASTA text and runs the ensemble
pipeline on every sequence. Results are written incrementally to TSV.

Usage (unchanged, fully backwards-compatible):
    from batch.processor import run_batch
    df = run_batch("proteome.fasta", evalue=1e-3, use_structure=True)

New: job management API
    from batch.processor import BatchManager
    mgr = BatchManager()
    job_id = mgr.submit("proteome.fasta", use_structure=True)
    print(mgr.status(job_id))   # queued | running | done | failed

Limits (all overridable via environment variables):
    ANDOM_MAX_SEQ   max sequences per batch job   (default 500)
    ANDOM_MAX_LEN   max amino acids per sequence  (default 5000)
    ANDOM_MIN_LEN   min amino acids per sequence  (default 20)
    ANDOM_MAX_JOBS  max concurrent jobs           (default 4)
"""
from __future__ import annotations

import os
import json
import uuid
import time
import threading
import pandas as pd
from pathlib import Path
from typing import Generator
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# ── configurable limits ───────────────────────────────────────────────────────
MAX_SEQUENCES: int = int(os.environ.get("ANDOM_MAX_SEQ",  500))
MAX_SEQ_LEN:   int = int(os.environ.get("ANDOM_MAX_LEN", 5000))
MIN_SEQ_LEN:   int = int(os.environ.get("ANDOM_MIN_LEN",   20))
MAX_JOBS:      int = int(os.environ.get("ANDOM_MAX_JOBS",    4))

# ── FASTA helpers ─────────────────────────────────────────────────────────────

def _parse_fasta(source: str) -> Generator:
    """
    Yield (seq_id, sequence) from a FASTA file path OR raw FASTA text.
    Accepts both multi-FASTA and bare sequences (no header).
    """
    if "\n" not in source and Path(source).exists():
        text = Path(source).read_text()
    else:
        text = source

    # bare sequence with no header
    lines = text.strip().splitlines()
    if lines and not lines[0].startswith(">"):
        text = f">query\n{text}"
        lines = text.splitlines()

    seq_id, buf = None, []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id:
                yield seq_id, "".join(buf)
            seq_id = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.upper())
    if seq_id:
        yield seq_id, "".join(buf)


def parse_fasta_text(text: str) -> dict[str, str]:
    """Parse raw FASTA text → {seq_id: sequence}. Used by the Streamlit UI."""
    return {sid: seq for sid, seq in _parse_fasta(text)}


def validate_sequences(seqs: dict[str, str]) -> list[str]:
    """
    Return a list of human-readable error strings.
    Empty list means all sequences are valid.
    """
    errors: list[str] = []
    valid_aa = set("ACDEFGHIKLMNPQRSTVWYBXZUOJ")

    if not seqs:
        errors.append("No sequences found in input.")
        return errors

    if len(seqs) > MAX_SEQUENCES:
        errors.append(
            f"{len(seqs)} sequences submitted — limit is {MAX_SEQUENCES}. "
            "Split your input or contact the administrator."
        )

    for sid, seq in seqs.items():
        if len(seq) < MIN_SEQ_LEN:
            errors.append(f"'{sid}': too short ({len(seq)} aa, minimum {MIN_SEQ_LEN}).")
        if len(seq) > MAX_SEQ_LEN:
            errors.append(f"'{sid}': too long ({len(seq)} aa, maximum {MAX_SEQ_LEN}).")
        bad = set(seq) - valid_aa
        if bad:
            errors.append(f"'{sid}': invalid characters {sorted(bad)}.")

    return errors


# ── original run_batch (unchanged public API) ─────────────────────────────────

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
    max_sequences   : stop after N sequences (None = all, hard cap MAX_SEQUENCES)

    Returns
    -------
    Combined DataFrame of all hits with query_id and source columns.
    """
    from search import sequence as seq_search
    from search import structure as str_search

    limit = min(max_sequences or MAX_SEQUENCES, MAX_SEQUENCES)
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    all_results: list[pd.DataFrame] = []
    header_written = False

    for i, (seq_id, seq) in enumerate(_parse_fasta(fasta_path)):
        if i >= limit:
            print(f"  Reached limit of {limit} sequences — stopping.")
            break

        if len(seq) < MIN_SEQ_LEN or len(seq) > MAX_SEQ_LEN:
            print(f"  [{i+1}] {seq_id} skipped (length {len(seq)} aa out of range)")
            continue

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


# ── Job management ────────────────────────────────────────────────────────────

class BatchManager:
    """
    Lightweight in-process job manager for the Streamlit batch UI.

    Jobs run in background threads. Status and results are kept in memory
    and also persisted to output/results/<job_id>/job.json so a page
    refresh doesn't lose them.

    Example
    -------
    mgr = BatchManager()
    job_id = mgr.submit("proteome.fasta", use_structure=True)
    print(mgr.status(job_id))        # "queued" → "running" → "done"
    df = mgr.results(job_id)         # pd.DataFrame once done
    mgr.cancel(job_id)               # request cancellation
    """

    def __init__(self, output_root: str = "output"):
        self._root    = Path(output_root) / "results"
        self._root.mkdir(parents=True, exist_ok=True)
        self._jobs: dict[str, dict] = {}          # job_id → state dict
        self._cancel: dict[str, bool] = {}        # job_id → cancel flag
        self._lock = threading.Lock()
        self._semaphore = threading.Semaphore(MAX_JOBS)
        self._reload_from_disk()

    # ── public API ─────────────────────────────────────────────────────────

    def submit(
        self,
        fasta_path: str,
        evalue: float = 1e-3,
        iterations: int = 3,
        use_structure: bool = False,
    ) -> str:
        """Enqueue a batch job. Returns job_id immediately."""
        job_id = uuid.uuid4().hex[:10]
        job_dir = self._root / job_id
        job_dir.mkdir(parents=True, exist_ok=True)

        state = {
            "job_id":      job_id,
            "fasta_path":  str(fasta_path),
            "status":      "queued",
            "submitted":   time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "started":     None,
            "finished":    None,
            "progress":    0,
            "total":       0,
            "error":       None,
            "out_tsv":     str(job_dir / "results.tsv"),
            "evalue":      evalue,
            "iterations":  iterations,
            "use_structure": use_structure,
        }

        with self._lock:
            self._jobs[job_id] = state
            self._cancel[job_id] = False
        self._save(state)

        t = threading.Thread(target=self._run, args=(job_id,), daemon=True)
        t.start()
        return job_id

    def status(self, job_id: str) -> str:
        """Return status string: queued | running | done | failed | cancelled."""
        with self._lock:
            return self._jobs.get(job_id, {}).get("status", "unknown")

    def info(self, job_id: str) -> dict:
        """Return full state dict for a job."""
        with self._lock:
            return dict(self._jobs.get(job_id, {}))

    def results(self, job_id: str) -> pd.DataFrame:
        """Return results DataFrame (empty if not done yet)."""
        info = self.info(job_id)
        tsv = info.get("out_tsv")
        if tsv and Path(tsv).exists():
            try:
                return pd.read_csv(tsv, sep="\t")
            except Exception:
                pass
        return pd.DataFrame()

    def cancel(self, job_id: str) -> None:
        """Request cancellation of a running or queued job."""
        with self._lock:
            self._cancel[job_id] = True
            if self._jobs.get(job_id, {}).get("status") == "queued":
                self._jobs[job_id]["status"] = "cancelled"
                self._save(self._jobs[job_id])

    def list_jobs(self) -> list[dict]:
        """All jobs, newest first."""
        with self._lock:
            jobs = list(self._jobs.values())
        return sorted(jobs, key=lambda j: j.get("submitted", ""), reverse=True)

    # ── internal ────────────────────────────────────────────────────────────

    def _run(self, job_id: str) -> None:
        from search import sequence as seq_search
        from search import structure as str_search

        with self._semaphore:
            with self._lock:
                state = self._jobs[job_id]
                if state["status"] == "cancelled":
                    return
                state["status"]  = "running"
                state["started"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            self._save(state)

            try:
                seqs = list(_parse_fasta(state["fasta_path"]))
                limit = min(len(seqs), MAX_SEQUENCES)
                all_results: list[pd.DataFrame] = []
                header_written = False
                out_tsv = state["out_tsv"]

                with self._lock:
                    state["total"] = limit

                for i, (seq_id, seq) in enumerate(seqs[:limit]):
                    if self._cancel.get(job_id):
                        with self._lock:
                            state["status"] = "cancelled"
                        self._save(state)
                        return

                    if not (MIN_SEQ_LEN <= len(seq) <= MAX_SEQ_LEN):
                        continue

                    fasta = f">{seq_id}\n{seq}"
                    df_seq, err = seq_search.run(
                        fasta,
                        evalue=state["evalue"],
                        iterations=state["iterations"],
                    )
                    if df_seq is not None and len(df_seq) > 0:
                        df_seq = df_seq.copy()
                        df_seq["query_id"] = seq_id
                        df_seq["source"]   = "SCOPe_sequence"
                        all_results.append(df_seq)
                        df_seq.to_csv(out_tsv, sep="\t", index=False,
                                      mode="a", header=not header_written)
                        header_written = True

                    if state["use_structure"] and len(seq) <= 400:
                        df_str, err2 = str_search.run(seq)
                        if df_str is not None and len(df_str) > 0:
                            df_str = df_str.copy()
                            df_str["query_id"] = seq_id
                            df_str["source"]   = "CATH_structure"
                            all_results.append(df_str)
                            df_str.to_csv(out_tsv, sep="\t", index=False,
                                          mode="a", header=not header_written)
                            header_written = True

                    with self._lock:
                        state["progress"] = i + 1
                    self._save(state)

                with self._lock:
                    state["status"]   = "done"
                    state["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                self._save(state)

            except Exception as exc:
                with self._lock:
                    state["status"] = "failed"
                    state["error"]  = str(exc)
                self._save(state)

    def _save(self, state: dict) -> None:
        path = self._root / state["job_id"] / "job.json"
        try:
            path.write_text(json.dumps(state, indent=2, default=str))
        except Exception:
            pass

    def _reload_from_disk(self) -> None:
        for job_dir in self._root.iterdir():
            job_file = job_dir / "job.json"
            if not job_file.exists():
                continue
            try:
                state = json.loads(job_file.read_text())
                if state.get("status") in ("running", "queued"):
                    state["status"] = "failed"
                    state["error"]  = "Server restarted — job did not complete."
                    job_file.write_text(json.dumps(state, indent=2))
                self._jobs[state["job_id"]] = state
                self._cancel[state["job_id"]] = False
            except Exception:
                pass
