FROM python:3.10-slim

LABEL maintainer="Dandekar Lab Wuerzburg"
LABEL description="AnDOM 2.0 — Structural Domain Finder"
LABEL version="2.0"

# ── system deps ───────────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl tar gzip \
    && rm -rf /var/lib/apt/lists/*

# ── working directory ─────────────────────────────────────────────────────────
WORKDIR /app

# ── python deps (cached layer) ────────────────────────────────────────────────
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ── install foldseek ──────────────────────────────────────────────────────────
RUN wget -q https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz \
    && tar xf foldseek-linux-avx2.tar.gz \
    && mv foldseek/bin/foldseek /usr/local/bin/foldseek \
    && rm -rf foldseek foldseek-linux-avx2.tar.gz

# ── install mmseqs2 ───────────────────────────────────────────────────────────
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
    && tar xf mmseqs-linux-avx2.tar.gz \
    && mv mmseqs/bin/mmseqs /usr/local/bin/mmseqs \
    && rm -rf mmseqs mmseqs-linux-avx2.tar.gz

# ── copy application code ─────────────────────────────────────────────────────
COPY app.py config.py ./
COPY src/ ./src/
COPY docs/ ./docs/

# ── create data and output mount points ───────────────────────────────────────
# data/ and output/ are NOT baked into the image — they are mounted at runtime
RUN mkdir -p /data /output/tmp /output/results

# ── streamlit config ──────────────────────────────────────────────────────────
RUN mkdir -p /root/.streamlit
COPY deploy/streamlit_config.toml /root/.streamlit/config.toml

# ── expose port ───────────────────────────────────────────────────────────────
EXPOSE 8501

# ── entrypoint ────────────────────────────────────────────────────────────────
CMD ["streamlit", "run", "app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true"]
