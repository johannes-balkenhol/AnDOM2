"""
Batch processing for AnDOM 2.0.

Multi-user safe: every job gets an isolated output directory.
MMseqs2/Foldseek intermediate files are written to output/tmp/<job_id>/
so the project root stays clean regardless of how many users run searches.

Retention / cleanup: run cleanup.sh via cron (daily recommended).
Environment variables:
    ANDOM_MAX_SEQ   max sequences per batch job   (default 500)
    ANDOM_MAX_LEN   max amino acids per sequence  (default 5000)
    ANDOM_MIN_LEN   min amino acids per sequence  (default 20)
    ANDOM_MAX_JOBS  max concurrent jobs           (default 4)
    ANDOM_OUTPUT_DIR  output root                 (default output/)
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

_OUTPUT_ROOT = os.environ.get("ANDOM_OUTPUT_DIR", "output")

# ── FASTA helpers ─────────────────────────────────────────────────────────────

def _parse_fasta(source: str) -> Generator:
    """Yield (seq_id, sequence) from a FASTA file path OR raw FASTA text."""
    if "\n" not in source and Path(source).exists():
        text = Path(source).read_text()
    else:
        text = source
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
    """Return list of error strings. Empty = all valid."""
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


# ── original run_batch (backwards-compatible public API) ──────────────────────

def run_batch(
    fasta_path: str,
    evalue: float = 1e-3,
    iterations: int = 3,
    use_structure: bool = False,
    out_tsv: str = "batch_results.tsv",
    max_sequences: int | None = None,
    job_tmp_dir: str | None = None,   # NEW: isolated tmp dir per job
) -> pd.DataFrame:
    """
    Run AnDOM 2.0 ensemble pipeline on every sequence in a FASTA file.
    All MMseqs2/Foldseek intermediate files go to job_tmp_dir (not project root).
    """
    from search import sequence as seq_search
    from search import structure as str_search

    limit = min(max_sequences or MAX_SEQUENCES, MAX_SEQUENCES)
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)

    # pass tmp dir to search modules so they don't pollute project root
    search_kwargs = {}
    if job_tmp_dir:
        search_kwargs["tmp_dir"] = job_tmp_dir

    all_results: list[pd.DataFrame] = []
    header_written = False

    for i, (seq_id, seq) in enumerate(_parse_fasta(fasta_path)):
        if i >= limit:
            print(f"  Reached limit of {limit} sequences — stopping.")
            break
        if not (MIN_SEQ_LEN <= len(seq) <= MAX_SEQ_LEN):
            print(f"  [{i+1}] {seq_id} skipped (length {len(seq)} aa)")
            continue

        print(f"  [{i+1}] {seq_id} ({len(seq)} aa)")
        fasta = f">{seq_id}\n{seq}"

        df_seq, err = seq_search.run(fasta, evalue=evalue, iterations=iterations,
                                     **search_kwargs)
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
            df_str, err2 = str_search.run(seq, **search_kwargs)
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

    return pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()


# ── BatchManager: multi-user job queue ────────────────────────────────────────

class BatchManager:
    """
    Thread-safe batch job manager.

    Multi-user design:
    - Each job gets output/results/<job_id>/ — isolated per job
    - MMseqs2 tmp files go to output/tmp/<job_id>/ — never in project root
    - Job state persisted to job.json — survives container restart
    - Completed/failed jobs older than JOB_DAYS cleaned by cleanup.sh
    """

    def __init__(self, output_root: str = _OUTPUT_ROOT):
        self._root    = Path(output_root) / "results"
        self._tmp     = Path(output_root) / "tmp"
        self._root.mkdir(parents=True, exist_ok=True)
        self._tmp.mkdir(parents=True, exist_ok=True)
        self._jobs: dict[str, dict] = {}
        self._cancel: dict[str, bool] = {}
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
        session_id: str | None = None,   # optional browser session tag
    ) -> str:
        job_id  = uuid.uuid4().hex[:10]
        job_dir = self._root / job_id
        tmp_dir = self._tmp  / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)

        state = {
            "job_id":       job_id,
            "session_id":   session_id or "anonymous",
            "fasta_path":   str(fasta_path),
            "status":       "queued",
            "submitted":    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "started":      None,
            "finished":     None,
            "progress":     0,
            "total":        0,
            "error":        None,
            "out_tsv":      str(job_dir / "results.tsv"),
            "job_dir":      str(job_dir),
            "tmp_dir":      str(tmp_dir),
            "evalue":       evalue,
            "iterations":   iterations,
            "use_structure":use_structure,
        }
        with self._lock:
            self._jobs[job_id] = state
            self._cancel[job_id] = False
        self._save(state)

        t = threading.Thread(target=self._run, args=(job_id,), daemon=True)
        t.start()
        return job_id

    def status(self, job_id: str) -> str:
        with self._lock:
            return self._jobs.get(job_id, {}).get("status", "unknown")

    def info(self, job_id: str) -> dict:
        with self._lock:
            return dict(self._jobs.get(job_id, {}))

    def results(self, job_id: str) -> pd.DataFrame:
        info = self.info(job_id)
        tsv = info.get("out_tsv")
        if tsv and Path(tsv).exists():
            try:
                return pd.read_csv(tsv, sep="\t")
            except Exception:
                pass
        return pd.DataFrame()

    def cancel(self, job_id: str) -> None:
        with self._lock:
            self._cancel[job_id] = True
            if self._jobs.get(job_id, {}).get("status") == "queued":
                self._jobs[job_id]["status"] = "cancelled"
                self._save(self._jobs[job_id])

    def list_jobs(self, session_id: str | None = None) -> list[dict]:
        """Return all jobs, or filtered by session_id if provided."""
        with self._lock:
            jobs = list(self._jobs.values())
        if session_id:
            jobs = [j for j in jobs if j.get("session_id") == session_id]
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
                seqs  = list(_parse_fasta(state["fasta_path"]))
                limit = min(len(seqs), MAX_SEQUENCES)
                all_results: list[pd.DataFrame] = []
                header_written = False
                out_tsv = state["out_tsv"]
                tmp_dir = state["tmp_dir"]

                with self._lock:
                    state["total"] = limit

                search_kwargs = {"tmp_dir": tmp_dir}

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
                        fasta, evalue=state["evalue"],
                        iterations=state["iterations"], **search_kwargs
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
                        df_str, err2 = str_search.run(seq, **search_kwargs)
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

                # clean up tmp dir for this job immediately after success
                import shutil
                try:
                    shutil.rmtree(state["tmp_dir"], ignore_errors=True)
                except Exception:
                    pass

            except Exception as exc:
                with self._lock:
                    state["status"] = "failed"
                    state["error"]  = str(exc)
                self._save(state)

    def _save(self, state: dict) -> None:
        path = Path(state["job_dir"]) / "job.json"
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
