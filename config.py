"""
config.py
=========
CarboxyDB — Central configuration. Every script imports from here.

Usage:
    from config import CFG, PATHS, LOG

Do NOT hardcode paths, thresholds, or API credentials in individual scripts.
Change them here once and all scripts pick up the new values.
"""

import logging
import os
from datetime import datetime
from pathlib import Path

# ── Project root (one level above scripts/) ───────────────────────────────────
# Works whether you run from the project root or from scripts/
_THIS_FILE = Path(__file__).resolve()
ROOT = _THIS_FILE.parent   # config.py lives at project root

# ── Timestamp (set once per import, shared by all scripts in a run) ───────────
TS = datetime.now().strftime("%Y%m%d_%H%M%S")


# ═══════════════════════════════════════════════════════════════════════════════
# PATHS  — all file/directory locations
# ═══════════════════════════════════════════════════════════════════════════════
class PATHS:
    # Raw downloads
    RAW             = ROOT / "data" / "raw"
    RAW_BRENDA      = ROOT / "data" / "raw" / "brenda"
    RAW_SWISSPROT   = ROOT / "data" / "raw" / "uniprot" / "swissprot"
    RAW_TREMBL      = ROOT / "data" / "raw" / "uniprot" / "trembl"
    RAW_NEGATIVES   = ROOT / "data" / "raw" / "uniprot" / "negatives"

    # Intermediate (per-step outputs, not final)
    INTERIM         = ROOT / "data" / "interim"

    # Primary merged dataset — single source of truth
    PRIMARY         = ROOT / "data" / "primary"
    MASTER_TSV      = ROOT / "data" / "primary" / "master.tsv"
    MASTER_FASTA    = ROOT / "data" / "primary" / "master.fasta"
    ID_MAP          = ROOT / "data" / "primary" / "id_map.tsv"   # CDB_ID ↔ uniprot_id

    # Feature layers (one TSV per layer, index = CDB_ID)
    FEATURES        = ROOT / "data" / "features"
    FEAT_COMP       = ROOT / "data" / "features" / "composition"
    FEAT_DOMAINS    = ROOT / "data" / "features" / "domains"
    FEAT_MOTIFS     = ROOT / "data" / "features" / "motifs"
    FEAT_BLAST      = ROOT / "data" / "features" / "blast"
    FEAT_ESM2       = ROOT / "data" / "features" / "esm2"
    FEAT_ANKH       = ROOT / "data" / "features" / "ankh"   # Ankh embeddings
    FEAT_INTERPRO   = ROOT / "data" / "features" / "interpro"  # InterProScan
    FEAT_MEME       = ROOT / "data" / "features" / "meme"
    MEME_HITS_TSV   = ROOT / "data" / "features" / "meme" / "meme_hits.tsv"  # PENDING

    # ML-ready splits
    ML              = ROOT / "data" / "ml"

    # Benchmark
    BENCHMARK       = ROOT / "data" / "benchmark"

    # SHAP outputs
    SHAP            = ROOT / "data" / "shap"

    # Reference databases
    PFAM_HMM        = ROOT / "data" / "dbs" / "pfam" / "Pfam-A.hmm"
    BLAST_DB        = ROOT / "data" / "dbs" / "blast" / "carboxy_train"
    PROSITE_DAT     = ROOT / "data" / "dbs" / "prosite" / "prosite.dat"

    # Models
    MODELS          = ROOT / "models"

    # SQLite database
    DB              = ROOT / "database" / "carbodb.sqlite"

    # Logs
    LOGS            = ROOT / "logs"


# ═══════════════════════════════════════════════════════════════════════════════
# CFG  — all tunable parameters
# ═══════════════════════════════════════════════════════════════════════════════
class CFG:

    # ── ID system ──────────────────────────────────────────────────────────────
    # CDB_ID format: CDB000001, CDB000002, …
    # Assigned once in script 03_merge_all_sources.py and never changed.
    CDB_ID_PREFIX   = "CDB"
    CDB_ID_WIDTH    = 6          # zero-padded to 6 digits → CDB000001

    # ── Sequence validation ────────────────────────────────────────────────────
    SEQ_MIN_LEN     = 50
    SEQ_MAX_LEN     = 5_000
    SEQ_VALID_AA    = set("ACDEFGHIKLMNPQRSTVWY")
    SEQ_MIN_UNIQUE  = 5          # reject sequences with < 5 distinct AA

    # ── UniProt download scale ─────────────────────────────────────────────────
    UNIPROT_BATCH   = 500        # rows per API page (max allowed)
    SP_MAX_PER_EC   = 100_000    # SwissProt: effectively unlimited per EC
    TB_MAX_PER_EC   = 100_000    # TrEMBL: cap per EC to stay manageable
    NEG_RATIO       = 5.0        # target negatives = NEG_RATIO × positives

    # ── BRENDA SOAP API ────────────────────────────────────────────────────────
    BRENDA_EMAIL    = os.environ.get("BRENDA_EMAIL", "")       # set in env
    BRENDA_PASSWORD = os.environ.get("BRENDA_PASSWORD", "")    # SHA256 hash
    BRENDA_WSDL     = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

    # ── CD-HIT clustering ──────────────────────────────────────────────────────
    CDHIT_THRESHOLD = 0.90       # 90% identity = same cluster
    CDHIT_WORDSIZE  = 5
    CDHIT_THREADS   = 8

    # ── ML data split ──────────────────────────────────────────────────────────
    SPLIT_TRAIN     = 0.70
    SPLIT_VAL       = 0.15
    SPLIT_TEST      = 0.15
    SPLIT_SEED      = 42

    # ── XGBoost hyper-parameters ───────────────────────────────────────────────
    XGB_PARAMS = dict(
        n_estimators          = 500,
        max_depth             = 6,
        learning_rate         = 0.05,
        subsample             = 0.8,
        colsample_bytree      = 0.8,
        min_child_weight      = 3,
        tree_method           = "hist",
        early_stopping_rounds = 30,
        verbosity             = 0,
        random_state          = 42,
    )

    # ── Feature layers ─────────────────────────────────────────────────────────
    # v3 = composition + domains + blast (no ESM-2)
    # v5 = v3 + ESM-2 embeddings
    FEATURE_VERSIONS = ["v3", "v5"]
    ESM2_MODEL       = "esm2_t33_650M_UR50D"
    ESM2_DIM         = 1280

    # ── HMMER / Pfam ───────────────────────────────────────────────────────────
    HMMER_THREADS    = 8
    HMMER_EVALUE     = 1e-3

    # ── BLAST ──────────────────────────────────────────────────────────────────
    BLAST_THREADS    = 8
    BLAST_EVALUE     = 0.001
    BLAST_MAX_HITS   = 1

    # ── Evidence tiers ─────────────────────────────────────────────────────────
    # Lower = higher confidence
    TIER_EXPERIMENTAL = 1   # BRENDA Km measured
    TIER_CURATED      = 2   # SwissProt reviewed
    TIER_PREDICTED    = 3   # TrEMBL / model output
    TIER_INFERRED     = 4   # BLAST / Pfam

    # ── Km ─────────────────────────────────────────────────────────────────────
    KM_UNIT          = "mM"
    KM_MIN_VALID     = 1e-5   # discard Km < 0.00001 mM as likely error
    KM_MAX_VALID     = 1_000  # discard Km > 1000 mM as likely error

    # ── CO2 substrate synonyms (for BRENDA parsing) ────────────────────────────
    CO2_SUBSTRATES = {
        "CO2", "co2", "carbon dioxide",
        "HCO3-", "HCO3", "hco3", "bicarbonate",
        "CO2/HCO3-", "CO(2)", "carbon-dioxide",
        "dissolved CO2",
    }

    # ── Confirmed carboxylase Pfam domains ─────────────────────────────────────
    CARBOXY_PFAM = {
        "PF00016", "PF02788", "PF00101", "PF00194", "PF03119",
        "PF00311", "PF00821", "PF02785", "PF00364", "PF01039",
        "PF02786", "PF02787", "PF00289", "PF01309", "PF03599",
        "PF03590", "PF00384", "PF00682",
    }

    # ── PROSITE patterns ───────────────────────────────────────────────────────
    PROSITE_IDS = [
        "PS00157", "PS00158", "PS00162", "PS00188",
        "PS00781", "PS00393", "PS00017",
        "PS00488", "PS00546", "PS01163",
        "PS00014", "PS00013", "PS00012", "PS00011",
    ]

    # ── Label system ───────────────────────────────────────────────────────────
    # Three-class labelling — scientifically honest, avoids FP contamination.
    # Script 03 assigns labels automatically based on EC number.
    #
    #   label=1  TRUE CARBOXYLASE    — fixes or interconverts CO2 as primary rxn
    #            used for: binary classification, EC prediction, Km regression
    #
    #   label=2  CO2-RELATED         — ancestral CO2-handling, not primary carbox.
    #            e.g. decarboxylases, dehydrogenases with CO2-binding motifs
    #            used for: excluded from binary training OR as hard negatives
    #            NOT used for Km regression
    #
    #   label=0  TRUE NEGATIVE       — no CO2 involvement at all
    #            used for: binary classification training
    #
    # Background: CarboxyPred v3 showed 70% of false positives were ICDH
    # (1.1.1.42) and MDH (1.1.1.37) — enzymes with ancestral CO2-binding motifs.
    # Forcing them into label=0 added noise. label=2 captures the biology.

    LABEL_POSITIVE   = 1   # true carboxylase
    LABEL_ANCESTRAL  = 2   # CO2-related but not primary carboxylase
    LABEL_NEGATIVE   = 0   # true negative

    # EC classes confirmed as ancestral/CO2-related (label=2)
    # These have CO2/HCO3- Km entries in BRENDA but are NOT true carboxylases.
    # Source: CarboxyPred v3 FP analysis (false_positive_breakdown, Dec 2025)
    # 70.1% ICDH, 18.7% MDH, 7.0% Malic enzyme, 3.5% 6-PGD, 0.5% other
    CO2_RELATED_EC = {
        # Decarboxylases / CO2-releasing enzymes
        "1.1.1.42",   # Isocitrate DH (ICDH) — 70% of FPs in v3
        "1.1.1.37",   # Malate DH (MDH)       — 18.7% of FPs
        "1.1.1.40",   # Malic enzyme           — 7.0% of FPs
        "1.1.1.44",   # 6-PGD                  — 3.5% of FPs
        "1.1.1.43",   # Phosphogluconate DH
        "1.2.4.1",    # Pyruvate DH (CO2 release)
        "1.2.4.2",    # Oxoglutarate DH
        "4.1.1.1",    # Pyruvate decarboxylase
        "4.1.1.3",    # Oxaloacetate decarboxylase
        "4.1.1.7",    # Benzoylformate decarboxylase
        "4.1.1.11",   # Aspartate 1-decarboxylase
        "4.1.1.15",   # Glutamate decarboxylase
        "4.1.1.17",   # Ornithine decarboxylase
        "4.1.1.18",   # Lysine decarboxylase
        "4.1.1.19",   # Arginine decarboxylase
        "4.1.1.20",   # Diaminopimelate decarboxylase
        "4.1.1.25",   # Tyrosine decarboxylase
        "4.1.1.28",   # DOPA decarboxylase
        "4.1.1.33",   # Diphosphomevalonate decarboxylase
        "4.1.1.36",   # PABA decarboxylase
        "4.1.1.48",   # Indole-3-pyruvate decarboxylase
        "4.1.1.50",   # Adenosylmethionine decarboxylase
        "1.17.1.9",   # Found in current script 01 run — verify after completion
        "1.17.1.10",  # Found in current script 01 run — verify after completion
    }

    # EC classes dynamically discovered by script 01 that need manual review
    # after each BRENDA run — add to CO2_EC_CONFIRMED or CO2_RELATED_EC
    CO2_EC_NEEDS_REVIEW = {
        "1.17.98.4", "1.17.98.5",  # seen in current run
        "2.5.1.143",                # seen in current run
        "1.2.7.7",                  # seen in current run
        "1.3.1.85",                 # seen in current run
        "1.4.4.2",                  # seen in current run
    }

    # ── HPC configuration ──────────────────────────────────────────────────────
    # Cluster: wbbi206 login node + SLURM compute nodes
    # GPU: NVIDIA A100 80GB PCIe on wbbi206 (NOT in SLURM — run directly)
    # CPU compute: wbbi201 (80 cores, 515GB), wbbi202 (56c, 257GB),
    #              wbbi203 (72c, 515GB, sometimes busy)
    # SLURM partitions: artemis (default), hades, apollo — all infinite timelimit
    # Best CPU node for annotation jobs: wbbi201 via hades partition

    HPC_GPU_NODE        = "wbbi206"       # login node with A100
    HPC_GPU_PARTITION   = None            # A100 NOT in SLURM — run in screen
    HPC_CPU_PARTITION   = "hades"         # best for CPU array jobs
    HPC_CPU_BEST_NODE   = "wbbi201"       # 80 cores, 515 GB
    HPC_MAX_CPU         = 80             # safe max for SLURM jobs on wbbi201
    HPC_A100_VRAM_GB    = 80

    # ESM-2 on A100 — run directly on wbbi206, not via SLURM
    ESM2_BATCH_SIZE     = 128            # safe for A100 80GB
    ESM2_RUN_ON_LOGIN   = True           # flag for script 04c

    # Ankh on CPU — SLURM array job on hades partition
    ANKH_MODEL          = "ankh-large"   # 1024-dim, CPU-runnable
    ANKH_DIM            = 1024
    ANKH_BATCH_SIZE     = 32             # per worker
    ANKH_THREADS        = 80            # use all of wbbi201

    # HMMER — SLURM array job
    HMMER_THREADS       = 80            # override base — use full node
    HMMER_EVALUE        = 1e-3

    # InterProScan — installed locally
    INTERPROSCAN_SH     = (ROOT / "data" / "dbs" / "interpro" /
                           "interproscan-5.72-103.0" / "interproscan.sh")
    INTERPROSCAN_THREADS = 80
    INTERPROSCAN_APPS   = (
        "Pfam,TIGRFAM,SUPERFAMILY,PANTHER,Gene3D,CDD,ProSiteProfiles"
    )

    # CD-HIT — single SLURM job
    CDHIT_THREADS       = 80            # override base

    # ── Database schema version ────────────────────────────────────────────────
    SCHEMA_VERSION   = "2.1"            # bumped for label=2 system


# ═══════════════════════════════════════════════════════════════════════════════
# Logging  — call setup_logging(name) at the top of each script
# ═══════════════════════════════════════════════════════════════════════════════

def setup_logging(script_name: str) -> logging.Logger:
    """
    Create a logger that writes to both stdout and a timestamped log file.

    Usage:
        from config import setup_logging
        log = setup_logging("03_merge")
    """
    PATHS.LOGS.mkdir(parents=True, exist_ok=True)
    log_file = PATHS.LOGS / f"{script_name}_{TS}.log"

    fmt = logging.Formatter(
        "%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logger = logging.getLogger(script_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # File handler
    fh = logging.FileHandler(log_file)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info("Log file: %s", log_file)
    return logger


# ═══════════════════════════════════════════════════════════════════════════════
# ID helpers
# ═══════════════════════════════════════════════════════════════════════════════

def make_cdb_id(n: int) -> str:
    """Format an integer as a CDB_ID string. make_cdb_id(1) → 'CDB000001'"""
    return f"{CFG.CDB_ID_PREFIX}{n:0{CFG.CDB_ID_WIDTH}d}"


def load_id_map() -> dict:
    """
    Load the CDB_ID ↔ uniprot_id mapping.
    Returns dict: {uniprot_id: cdb_id}
    Raises FileNotFoundError if id_map.tsv does not yet exist.
    """
    import pandas as pd
    df = pd.read_csv(PATHS.ID_MAP, sep="\t", dtype=str)
    return dict(zip(df["uniprot_id"], df["cdb_id"]))


def latest_file(directory: Path, pattern: str) -> Path:
    """
    Return the most recent file matching a glob pattern.
    Raises FileNotFoundError if nothing matches.

    Example:
        ec_file = latest_file(PATHS.RAW_BRENDA, "co2_ec_classes_*.txt")
    """
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"No file matching '{pattern}' in {directory}"
        )
    return files[-1]

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
