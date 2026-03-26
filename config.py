"""
AnDOM 2.0 — central configuration.
All paths, constants, and example sequences live here.
Import this module from anywhere in the project.
"""
from pathlib import Path

# ── project root & key directories ────────────────────────────────────────────
import os
ROOT_DIR   = Path(__file__).parent
# allow Docker to override data/output paths via environment variables
_DATA_OVERRIDE   = os.environ.get("ANDOM_DATA_DIR")
_OUTPUT_OVERRIDE = os.environ.get("ANDOM_OUTPUT_DIR")
SRC_DIR    = ROOT_DIR / "src"
DATA_DIR   = Path(_DATA_OVERRIDE) if _DATA_OVERRIDE else ROOT_DIR / "data"
DOCS_DIR   = ROOT_DIR / "docs"
OUTPUT_DIR = Path(_OUTPUT_OVERRIDE) if _OUTPUT_OVERRIDE else ROOT_DIR / "output"
TMP_DIR    = str(OUTPUT_DIR / "tmp")

# ── external tool paths ───────────────────────────────────────────────────────
import shutil as _shutil
_mmseqs_local = DATA_DIR / "mmseqs" / "bin" / "mmseqs"
MMSEQS = str(_mmseqs_local) if _mmseqs_local.exists() else (_shutil.which("mmseqs") or "mmseqs")
FOLDSEEK = "foldseek"

# ── database paths ────────────────────────────────────────────────────────────
SCOPE_DB = str(DATA_DIR / "scopeSeqDB")
SCOPE_FA = str(DATA_DIR / "scopeseq_40.fa")
CATH_DB  = str(DATA_DIR / "cathDB")

# ── external APIs ─────────────────────────────────────────────────────────────
ESMFOLD_API    = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAXLEN = 400

# ── SCOP visual constants ─────────────────────────────────────────────────────
SCOP_COLORS = {
    "a": "#e74c3c", "b": "#3498db", "c": "#2ecc71",
    "d": "#f39c12", "e": "#9b59b6", "f": "#1abc9c", "g": "#e67e22",
}
SCOP_CLASSES = {
    "a": "All alpha",   "b": "All beta",     "c": "Alpha/beta",
    "d": "Alpha+beta",  "e": "Multi-domain", "f": "Membrane",
    "g": "Small proteins",
}

# ── example sequences ─────────────────────────────────────────────────────────
EXAMPLES = {
    "Ex 1": {
        "label": "Ex 1 — Haemoglobin (globin, all-alpha)",
        "seq":   "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "desc":  "Human haemoglobin alpha (142 aa) — classic all-alpha globin fold.",
    },
    "Ex 2": {
        "label": "Ex 2 — HIV-1 protease (aspartyl, all-beta)",
        "seq":   "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
        "desc":  "HIV-1 protease (99 aa) — all-beta aspartyl protease, major drug target.",
    },
    "Ex 3": {
        "label": "Ex 3 — TIM barrel (alpha/beta classic)",
        "seq":   "MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLSQELGASNEILLGAQNVDLNLPKDKFVVLIIVYNKPKDILFKDNEALENGKPFKQNLKAKDLALAGVTPDKMKDLKAKGISGAFVPNIVNLHSQAPADCLMSKLVAGEFEGNIYMGLKPNPEELAAAKSSKLSELIQAAYATGNQVAFKPLTDAAQKAAQESSGKKSATIFAGQATVEDGDTVYL",
        "desc":  "Triosephosphate isomerase (248 aa) — canonical TIM barrel. Sequence-only search.",
    },
    "Ex 4": {
        "label": "Ex 4 — Src SH2+SH3 (multi-domain)",
        "seq":   "MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETSGY",
        "desc":  "Src kinase SH2+SH3 domains (180 aa) — multi-domain signalling protein.",
    },
    "Ex 5": {
        "label": "Ex 5 — p53 DNA-binding domain (beta-sandwich)",
        "seq":   "SVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL",
        "desc":  "p53 tumour suppressor DNA-binding domain (130 aa) — beta-sandwich.",
    },
}
