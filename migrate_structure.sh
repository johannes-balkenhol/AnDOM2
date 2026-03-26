#!/bin/bash
# AnDOM 2.0 — migrate to proper project folder structure
# Run from ~/Projects_shared/Andom
set -e
echo "Migrating to proper project structure..."

# ── 1. create new folders ──────────────────────────────────────────────────────
mkdir -p src/search src/db src/batch src/benchmark
mkdir -p data docs tests output/tmp output/results

# ── 2. move source modules into src/ ──────────────────────────────────────────
mv search/sequence.py  src/search/sequence.py
mv search/structure.py src/search/structure.py
mv search/ensemble.py  src/search/ensemble.py
mv search/__init__.py  src/search/__init__.py

mv db/lookup.py        src/db/lookup.py
mv db/__init__.py      src/db/__init__.py

mv batch/processor.py  src/batch/processor.py
mv batch/__init__.py   src/batch/__init__.py

mv benchmark/run.py    src/benchmark/run.py
mv benchmark/__init__.py src/benchmark/__init__.py

# copy test_sets and results structure
cp -r benchmark/test_sets src/benchmark/test_sets 2>/dev/null || mkdir -p src/benchmark/test_sets
cp -r benchmark/results   src/benchmark/results   2>/dev/null || mkdir -p src/benchmark/results

# remove old empty dirs
rm -rf search/ db/ batch/ benchmark/

# ── 3. move data files into data/ ─────────────────────────────────────────────
# SCOPe FASTA
mv scopeseq_40.fa data/ 2>/dev/null || true

# SCOPe MMseqs2 DB
for f in scopeSeqDB scopeSeqDB.dbtype scopeSeqDB_h scopeSeqDB_h.dbtype \
          scopeSeqDB_h.index scopeSeqDB.index scopeSeqDB.lookup scopeSeqDB.source; do
    mv "$f" data/ 2>/dev/null || true
done

# CATH Foldseek DB
for f in cathDB cathDB.dbtype cathDB.index cathDB.lookup cathDB.source cathDB.version \
          cathDB_ca cathDB_ca.dbtype cathDB_ca.index \
          cathDB_clu cathDB_clu.dbtype cathDB_clu.index \
          cathDB_h cathDB_h.dbtype cathDB_h.index \
          cathDB_mapping cathDB_ss cathDB_ss.dbtype cathDB_ss.index \
          cathDB_taxonomy \
          cathDB_seq.0 cathDB_seq.1 cathDB_seq.dbtype cathDB_seq.index \
          cathDB_seq.lookup cathDB_seq.source cathDB_seq_mapping cathDB_seq_taxonomy \
          cathDB_seq_ca.0 cathDB_seq_ca.1 cathDB_seq_ca.dbtype cathDB_seq_ca.index \
          cathDB_seq_h.0 cathDB_seq_h.1 cathDB_seq_h.dbtype cathDB_seq_h.index \
          cathDB_seq_ss.0 cathDB_seq_ss.1 cathDB_seq_ss.dbtype cathDB_seq_ss.index; do
    mv "$f" data/ 2>/dev/null || true
done

# mmseqs binary
mv mmseqs/ data/ 2>/dev/null || true
mv mmseqs-linux-avx2.tar.gz data/ 2>/dev/null || true

# move runtime output files
mv tmp/ output/ 2>/dev/null || true
mv *.tsv output/results/ 2>/dev/null || true
mv *.pdb output/ 2>/dev/null || true
mv queryDB* output/ 2>/dev/null || true
mv resultDB* output/ 2>/dev/null || true
mv query_input.fasta output/ 2>/dev/null || true

# ── 4. update config.py with new paths ────────────────────────────────────────
cat > config.py << 'EOF'
"""
AnDOM 2.0 — central configuration.
All paths, constants, and example sequences live here.
Import this module from anywhere in the project.
"""
from pathlib import Path

# ── project root & key directories ────────────────────────────────────────────
ROOT_DIR   = Path(__file__).parent
SRC_DIR    = ROOT_DIR / "src"
DATA_DIR   = ROOT_DIR / "data"
DOCS_DIR   = ROOT_DIR / "docs"
OUTPUT_DIR = ROOT_DIR / "output"
TMP_DIR    = str(OUTPUT_DIR / "tmp")

# ── external tool paths ───────────────────────────────────────────────────────
MMSEQS   = str(DATA_DIR / "mmseqs" / "bin" / "mmseqs")
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
EOF

# ── 5. update sys.path in all src modules ─────────────────────────────────────
# Replace the old sys.path.insert with the correct one pointing up two levels
for f in src/search/sequence.py src/search/structure.py src/search/ensemble.py \
         src/db/lookup.py src/batch/processor.py src/benchmark/run.py; do
    sed -i "s|sys.path.insert(0, str(Path(__file__).parent.parent))|sys.path.insert(0, str(Path(__file__).parent.parent.parent))|g" "$f"
done

# ── 6. update app.py to use src/ ──────────────────────────────────────────────
sed -i "s|sys.path.insert(0, str(Path(__file__).parent))|import sys\nfrom pathlib import Path\nsys.path.insert(0, str(Path(__file__).parent / 'src'))\nsys.path.insert(0, str(Path(__file__).parent))|g" app.py 2>/dev/null || true

# ── 7. create tests/ scaffold ─────────────────────────────────────────────────
touch tests/__init__.py

cat > tests/test_sequence.py << 'EOF'
"""Unit tests for search/sequence.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))

HBA_SEQ = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"

def test_fasta_format():
    fasta = f">HBA_HUMAN\n{HBA_SEQ}"
    assert fasta.startswith(">")
    assert len(HBA_SEQ) == 141

def test_scop_lookup_loads():
    import db.lookup as lk
    domains = lk.all_domains()
    assert len(domains) > 10000, "SCOPe lookup should have >10k domains"

def test_pdb_url():
    import db.lookup as lk
    url = lk.pdb_url("d3d1ka_")
    assert "3d1k" in url

if __name__ == "__main__":
    test_fasta_format()
    print("test_fasta_format passed")
    test_scop_lookup_loads()
    print("test_scop_lookup_loads passed")
    test_pdb_url()
    print("test_pdb_url passed")
    print("All tests passed.")
EOF

cat > tests/test_ensemble.py << 'EOF'
"""Unit tests for search/ensemble.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))
import pandas as pd

def test_domain_bar_empty():
    from search.ensemble import domain_bar_html
    html = domain_bar_html(None, None, 100)
    assert "<div" in html

def test_evalue_to_score():
    from search.ensemble import _evalue_to_score
    assert _evalue_to_score(1e-30) == 1.0
    assert _evalue_to_score(1.0) == 0.0
    assert 0 < _evalue_to_score(1e-10) < 1.0

if __name__ == "__main__":
    test_domain_bar_empty()
    print("test_domain_bar_empty passed")
    test_evalue_to_score()
    print("test_evalue_to_score passed")
    print("All tests passed.")
EOF

# ── 8. update .gitignore ───────────────────────────────────────────────────────
cat > .gitignore << 'EOF'
# ── data (too large for git, download separately) ─────────────────────────────
data/

# ── runtime output ─────────────────────────────────────────────────────────────
output/

# ── python ────────────────────────────────────────────────────────────────────
__pycache__/
*.pyc
*.pyo
.venv/
*.egg-info/

# ── streamlit ─────────────────────────────────────────────────────────────────
.streamlit/

# ── editor ────────────────────────────────────────────────────────────────────
.vscode/
.idea/
*.swp
EOF

# ── 9. update README ───────────────────────────────────────────────────────────
cat > README.md << 'EOF'
# AnDOM 2.0 — Structural Domain Finder

Protein domain assignment by ensemble sequence + structure search.
Update of Schmidt, Bork & Dandekar (*J Chem Inf Comput Sci* 2002).

## Project structure

```
AnDOM2/
├── app.py              # Streamlit entry point
├── config.py           # all paths, constants, examples
├── requirements.txt
├── src/
│   ├── search/         # sequence.py · structure.py · ensemble.py
│   ├── db/             # lookup.py (SCOPe metadata)
│   ├── batch/          # processor.py (multi-FASTA batch)
│   └── benchmark/      # run.py (coverage + precision benchmarks)
├── docs/               # SVG diagrams, method figures
├── tests/              # unit tests
├── data/               # databases (gitignored — download separately)
└── output/             # runtime temp files (gitignored)
```

## Install

```bash
mamba create -n andom python=3.10 -y
mamba activate andom
mamba install -c bioconda foldseek mmseqs2 -y
pip install -r requirements.txt
```

## Database setup

```bash
# SCOPe 2.08 ASTRAL sequences
wget --no-check-certificate \
  https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa \
  -O data/scopeseq_40.fa
mmseqs createdb data/scopeseq_40.fa data/scopeSeqDB

# CATH50 structural database
foldseek databases CATH50 data/cathDB output/tmp
```

## Run

```bash
streamlit run app.py --server.port 8501
```

## Batch processing

```python
from src.batch.processor import run_batch
df = run_batch("my_proteins.fasta", use_structure=True)
```

## Benchmark

```bash
python3 src/benchmark/run.py data/my_test_set.fasta --gold data/gold_standard.tsv
```

## References

- Schmidt S, Bork P, Dandekar T (2002) *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) *Science* 379:1123–1130
- Fox NK et al. (2014) *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) *Nucleic Acids Res* 49:D266–273
EOF

echo ""
echo "Migration complete. New structure:"
find . -name "*.py" -not -path "*/__pycache__/*" -not -path "./data/*" -not -path "./output/*" | sort
echo ""
echo "Next:"
echo "  1. python3 -c \"import sys; sys.path.insert(0,'src'); sys.path.insert(0,'.'); from config import EXAMPLES; import db.lookup as lk; print('OK', len(lk.all_domains()))\""
echo "  2. streamlit run app.py --server.port 8501"
echo "  3. git add -A && git commit -m 'refactor: proper project structure src/ data/ docs/ tests/ output/'"

# ── 10. rewrite app.py imports cleanly ────────────────────────────────────────
python3 << 'PYEOF'
content = open("app.py").read()

# Fix the sys.path to point at src/ for module imports
old_path = "import sys\nfrom pathlib import Path\nsys.path.insert(0, str(Path(__file__).parent / 'src'))\nsys.path.insert(0, str(Path(__file__).parent))"
new_path = """import sys
from pathlib import Path
_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "src"))"""

if "sys.path.insert(0, str(_ROOT / " not in content:
    content = content.replace(
        'import sys\nfrom pathlib import Path\nsys.path.insert(0, str(Path(__file__).parent))',
        new_path
    )
    open("app.py", "w").write(content)
    print("app.py imports updated")
else:
    print("app.py imports already correct")
PYEOF
