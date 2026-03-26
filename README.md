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
