#!/bin/bash
# Run once on the server to download databases into data/
# Usage: bash deploy/setup_data.sh
set -e

mkdir -p data output/tmp output/results

echo "=== Downloading SCOPe 2.08 ASTRAL sequences ==="
wget --no-check-certificate \
  "https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa" \
  -O data/scopeseq_40.fa

echo "=== Building MMseqs2 SCOPe database ==="
mmseqs createdb data/scopeseq_40.fa data/scopeSeqDB

echo "=== Downloading CATH50 Foldseek database ==="
foldseek databases CATH50 data/cathDB output/tmp --threads 8

echo ""
echo "Done. data/ is ready. Start the app with:"
echo "  docker compose up -d       # production"
echo "  streamlit run app.py       # development"
