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
