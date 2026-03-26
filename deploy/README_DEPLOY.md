# AnDOM 2.0 — Deployment Guide

## Option A: Docker (recommended for servers)

### 1. Clone the repo
```bash
git clone git@github.com:johannes-balkenhol/AnDOM2.git
cd AnDOM2
```

### 2. Download databases (once, ~1.5 GB)
```bash
bash deploy/setup_data.sh
```

### 3. Build and start
```bash
docker compose up -d
```

App is now at http://your-server:8501

### 4. View logs
```bash
docker compose logs -f andom
```

### 5. Update to new version
```bash
git pull
docker compose up -d --build
```

---

## Option B: Direct (development / local)

```bash
mamba create -n andom python=3.10 -y
mamba activate andom
mamba install -c bioconda foldseek mmseqs2 -y
pip install -r requirements.txt
bash deploy/setup_data.sh
streamlit run app.py --server.port 8501
```

---

## Option C: With nginx reverse proxy (HTTPS)

1. Install nginx and certbot
2. Copy `deploy/nginx.conf` to `/etc/nginx/sites-available/andom`
3. Edit the server_name to your domain
4. Run `certbot --nginx -d your-domain.de`
5. `docker compose up -d`

---

## Folder structure expected on server

```
AnDOM2/
├── app.py
├── config.py
├── docker-compose.yml
├── Dockerfile
├── src/
├── docs/
├── data/          ← created by setup_data.sh, gitignored
│   ├── scopeseq_40.fa
│   ├── scopeSeqDB*
│   └── cathDB*
└── output/        ← created at runtime, gitignored
    ├── tmp/
    └── results/
```

## Environment variables (optional override in docker-compose.yml)

| Variable | Default | Description |
|---|---|---|
| ANDOM_DATA_DIR | /data | Path to databases inside container |
| ANDOM_OUTPUT_DIR | /output | Path to runtime output inside container |
