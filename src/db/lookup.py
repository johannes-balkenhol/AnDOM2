"""
SCOPe and CATH domain metadata lookup.

Public API:
    get(domain_id)           -> SCOPe metadata dict
    all_domains()            -> full SCOPe lookup dict
    pdb_url(domain_id)       -> RCSB URL from SCOPe domain ID
    cath_pdb_url(target)     -> RCSB URL from CATH domain ID
    get_cath_code(domain)    -> CATH code from CATH domain ID e.g. "1.10.490.10"
    pdb_to_cath_code(pdb4)   -> CATH code from 4-char PDB code
    pdb_to_sccs(pdb4)        -> SCOPe sccs from 4-char PDB code
    pdb_to_uniprot(pdb4)     -> UniProt accession from PDB via SIFTS
    load_cath_codes()        -> full {domain_id: cath_code} dict
    cath_code_to_scop_class  -> CATH class digit -> SCOPe class letter
"""
import functools
import gzip
import streamlit as st
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import SCOPE_FA, SCOPE_FA_95, DATA_DIR


# ── SCOPe FASTA loading ───────────────────────────────────────────────────────

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


def get(domain_id: str) -> dict:
    """Return metadata dict for a SCOPe domain ID."""
    return _load().get(domain_id, {})


def all_domains() -> dict:
    """Return the full SCOPe lookup dict (40% + 95%)."""
    d = _load(SCOPE_FA)
    d.update(_load(SCOPE_FA_95))
    return d


# ── CATH domain list loading ──────────────────────────────────────────────────

@st.cache_data
def _load_cath(cath_file: str | None = None) -> dict:
    """
    Parse CathDomainList into {domain_id: cath_code}.
    Format: 1oaiA00  1  10  8  10  ...  → code = "1.10.8.10"
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
            if len(parts) < 5:
                continue
            domain_id = parts[0].lower()
            try:
                cath_code = f"{parts[1]}.{parts[2]}.{parts[3]}.{parts[4]}"
                result[domain_id] = cath_code
            except (IndexError, ValueError):
                pass
    return result


# ── PDB-code → CATH code ──────────────────────────────────────────────────────

@st.cache_data
def _pdb_to_cath_map() -> dict:
    """
    Build {pdb4: cath_code} from cath-domain-list.txt.
    First (lowest line number = best scoring) entry per PDB code wins.
    PDB code = first 4 chars of domain ID: '1oaiA00' -> '1oai'.
    """
    path = Path(DATA_DIR) / "cath-domain-list.txt"
    result: dict = {}
    if not path.exists():
        return result
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            pdb4 = parts[0][:4].lower()
            if pdb4 in result:
                continue
            try:
                result[pdb4] = f"{parts[1]}.{parts[2]}.{parts[3]}.{parts[4]}"
            except (IndexError, ValueError):
                pass
    return result


def pdb_to_cath_code(pdb4: str) -> str:
    """
    Return CATH code for a 4-character PDB code.
    e.g. '1c7d' -> '1.10.490.10'
    Returns '' if not found.
    """
    if not pdb4 or len(pdb4) < 4:
        return ""
    return _pdb_to_cath_map().get(pdb4.lower()[:4], "")


# ── PDB-code → SCOPe sccs ─────────────────────────────────────────────────────

@st.cache_data
def _pdb_to_sccs_map() -> dict:
    """
    Build {pdb4: sccs} from the SCOPe ASTRAL lookup.
    First SCOPe entry per PDB code wins.
    SCOPe domain ID format: d1hbaa_ -> PDB = chars [1:5] = '1hba'.
    """
    lk = all_domains()
    result: dict = {}
    for domain_id, info in lk.items():
        if len(domain_id) < 5:
            continue
        pdb4 = domain_id[1:5].lower()
        if pdb4 and pdb4 not in result:
            sccs = info.get("sccs", "")
            if sccs and sccs not in ("?", ""):
                result[pdb4] = sccs
    return result


def pdb_to_sccs(pdb4: str) -> str:
    """
    Return SCOPe sccs for a 4-character PDB code.
    e.g. '2dfp' -> 'a.1.1.2'
    Returns '' if not found in SCOPe ASTRAL 95%.
    """
    if not pdb4 or len(pdb4) < 4:
        return ""
    return _pdb_to_sccs_map().get(pdb4.lower()[:4], "")


# ── SIFTS: PDB → UniProt ──────────────────────────────────────────────────────

@st.cache_data
def _load_sifts() -> dict:
    """
    Load SIFTS mapping: {pdb4: uniprot_id} and {pdb4_chain: uniprot_id}.
    File: sifts_scop_cath.tsv.gz
    Columns: PDB  CHAIN  SF_DOMID  SP_PRIMARY  RES_BEG  RES_END  ...
    """
    sifts_file = Path(DATA_DIR) / "sifts_scop_cath.tsv.gz"
    result: dict = {}
    if not sifts_file.exists():
        return result
    with gzip.open(sifts_file, "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            pdb4  = parts[0].lower()
            chain = parts[1].upper()
            uid   = parts[3].strip()
            if not uid or uid == "-":
                continue
            key_chain = f"{pdb4}_{chain}"
            if key_chain not in result:
                result[key_chain] = uid
            if pdb4 not in result:
                result[pdb4] = uid
    return result


def pdb_to_uniprot(pdb4: str, chain: str = "") -> str:
    """
    Return UniProt accession for a PDB entry via SIFTS.
    e.g. pdb_to_uniprot('1hba', 'A') -> 'P68871'
    Falls back to pdb4-only if chain not found.
    Returns '' if not found.
    """
    sifts = _load_sifts()
    if chain:
        uid = sifts.get(f"{pdb4.lower()}_{chain.upper()}", "")
        if uid:
            return uid
    return sifts.get(pdb4.lower()[:4], "")


# ── public helpers ────────────────────────────────────────────────────────────

def get_cath_code(domain_id: str) -> str:
    """
    Return CATH code for a domain ID (CATH or SCOPe format).
    - CATH domain ID '1c7dA02'     -> direct lookup -> '1.10.490.10'
    - AlphaFold ID   'af_P9_1_136_1.10.490.10' -> embedded code
    - SCOPe domain   'd2dfoa_'     -> PDB code lookup -> first CATH entry
    - 4-char PDB     '1c7d'        -> first CATH entry for that PDB
    Returns '' if not found.
    """
    if not domain_id:
        return ""

    # Direct CATH domain ID lookup (7-char format: 1c7dA02)
    result = _load_cath().get(domain_id.lower(), "")
    if result:
        return result

    # AlphaFold entries embed CATH code: af_P9WN25_1_136_1.10.490.10
    if str(domain_id).startswith("af_"):
        parts = domain_id.split("_")
        for part in reversed(parts):
            if part.count(".") == 3:
                try:
                    [int(x) for x in part.split(".")]
                    return part
                except ValueError:
                    pass
        return ""

    # SCOPe domain ID: d2dfoa_ -> PDB = chars [1:5]
    if len(domain_id) >= 5 and domain_id[0] == "d":
        pdb4 = domain_id[1:5].lower()
        return pdb_to_cath_code(pdb4)

    # Plain 4-char PDB code
    if len(domain_id) == 4:
        return pdb_to_cath_code(domain_id)

    # Fallback: try first 4 chars
    return pdb_to_cath_code(domain_id[:4])


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
    Map CATH class digit to SCOPe class letter.
      1 = mainly alpha  -> a
      2 = mainly beta   -> b
      3 = alpha/beta    -> c
      4 = few secondary -> g
    """
    try:
        cath_class = cath_code.split(".")[0]
        return {"1": "a", "2": "b", "3": "c", "4": "g"}.get(cath_class, "")
    except Exception:
        return ""
