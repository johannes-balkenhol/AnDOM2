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
