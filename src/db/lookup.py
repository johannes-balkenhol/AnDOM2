"""
SCOPe and CATH domain metadata lookup.

Public API:
    get(domain_id)        -> SCOPe metadata dict
    all_domains()         -> full SCOPe lookup dict
    pdb_url(domain_id)    -> RCSB URL from SCOPe domain ID
    cath_pdb_url(target)  -> RCSB URL from CATH domain ID
    get_cath_code(domain) -> CATH classification code e.g. "1.10.490.10"
    load_cath_codes()     -> full {domain_id: cath_code} dict
"""
import streamlit as st
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import SCOPE_FA, SCOPE_FA_95, DATA_DIR


@st.cache_data
def _load(fa_file: str = SCOPE_FA, include_95: bool = True) -> dict:
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


@st.cache_data
def _load_cath(cath_file: str | None = None) -> dict:
    """
    Parse CathDomainList file into {domain_id: cath_code} dict.
    Format: 1oaiA00  1  10  8  10  ...
    CATH code = col0.col1.col2.col3 = "1.10.8.10"
    """
    if cath_file is None:
        cath_file = str(Path(DATA_DIR) / "cath-domain-list.txt")
    path = Path(cath_file)
    if not path.exists():
        return {}
    result: dict = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            domain_id = parts[0].lower()
            try:
                cath_code = f"{parts[1]}.{parts[2]}.{parts[3]}.{parts[4]}"
            except IndexError:
                cath_code = ".".join(parts[1:5])
            result[domain_id] = cath_code
    return result


def get(domain_id: str) -> dict:
    """Return metadata dict for a SCOPe domain ID."""
    return _load().get(domain_id, {})


def all_domains() -> dict:
    """Return the full SCOPe lookup dict (40% + 95%)."""
    d = _load(SCOPE_FA)
    d.update(_load(SCOPE_FA_95))
    return d


def get_cath_code(domain_id: str) -> str:
    """
    Return CATH classification code for a domain ID.
    e.g. '1c7dA02' -> '1.10.490.10'
    Domain ID is lowercased for lookup.
    Returns '' if not found.
    """
    # Try direct lookup first
    result = _load_cath().get(domain_id.lower(), "")
    if result:
        return result
    # AlphaFold entries embed CATH code in name: af_P9WN25_1_136_1.10.490.10
    if domain_id.startswith("af_"):
        parts = domain_id.split("_")
        for part in reversed(parts):
            if part.count(".") == 3:
                try:
                    [int(x) for x in part.split(".")]
                    return part
                except ValueError:
                    pass
    return ""


def load_cath_codes() -> dict:
    """Return full {domain_id: cath_code} dict."""
    return _load_cath()


def pdb_url(domain_id: str) -> str:
    """RCSB link from SCOPe domain ID (e.g. d3d1ka_ -> 3d1k)."""
    try:
        return f"https://www.rcsb.org/structure/{domain_id[1:5].lower()}"
    except Exception:
        return "https://scop.berkeley.edu"


def cath_pdb_url(target: str) -> str:
    """RCSB link from CATH domain ID (e.g. 1c7dA02 -> 1c7d)."""
    try:
        return f"https://www.rcsb.org/structure/{target[:4].lower()}"
    except Exception:
        return "https://www.cathdb.info"

def cath_code_to_scop_class(cath_code: str) -> str:
    """
    Map CATH classification code to approximate SCOPe class letter.
    Uses CATH class number (first digit of code):
      1 = mainly alpha    → a
      2 = mainly beta     → b
      3 = alpha/beta      → c  (mixed; SCOP separates c=barrel, d=orthogonal)
      4 = few secondary   → g  (small proteins)
    Returns '' if mapping not possible.
    """
    try:
        cath_class = cath_code.split(".")[0]
        return {"1": "a", "2": "b", "3": "c", "4": "g"}.get(cath_class, "")
    except Exception:
        return ""
