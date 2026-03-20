# AnDOM 2.0 — Structural Domain Finder

Protein domain assignment by ensemble sequence+structure search — SCOPe 2.08, CATH50, MMseqs2, ESMFold, Foldseek. Update of AnDOM 2002.

## What it does
Assigns structural domains to a protein sequence using two complementary layers:
- **Sequence search:** MMseqs2 PSI-search against SCOPe 2.08 ASTRAL (15,177 domains, 40% identity cutoff)
- **Structural search:** ESMFold structure prediction + Foldseek search against CATH50

Results are displayed as a two-row colour-coded domain architecture map with metadata, organism, fold classification, and direct PDB links.

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
wget --no-check-certificate https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa -O scopeseq_40.fa
mmseqs createdb scopeseq_40.fa scopeSeqDB

# CATH50 structural database
foldseek databases CATH50 cathDB tmp
```

## Run
```bash
streamlit run app.py --server.port 8501
```

## Reference
- Schmidt S, Bork P, Dandekar T (2002) A versatile structural domain analysis server using profile weight matrices. *J Chem Inf Comput Sci* 42:405-407
- AnDOM 2.0 — Dandekar Lab, University of Wuerzburg, 2026
