#!/bin/bash
# AnDOM 2.0 deploy script
# Usage:
#   ./deploy.sh prod    — pull main, rebuild prod container, restart prod only
#   ./deploy.sh dev     — start/restart dev container (uses live code mount, no rebuild)
#   ./deploy.sh dev-build — rebuild dev image (only needed after requirements change)
#   ./deploy.sh status  — show both containers
#   ./deploy.sh logs [prod|dev] — tail logs

set -e
cd /var/www/AnDOM2

case "$1" in

  prod)
    echo "=== Deploying PRODUCTION ==="
    echo "Pulling main branch..."
    git fetch origin main
    git checkout main
    git pull origin main
    echo "Rebuilding prod image..."
    docker compose -f docker-compose.yml build --no-cache
    echo "Restarting prod container (port 8501)..."
    docker compose -f docker-compose.yml down && docker compose -f docker-compose.yml up -d
    echo "✓ Production live at https://andom2.bioinfo-wuerz.eu"
    docker compose -f docker-compose.yml logs --tail=5
    ;;

  dev)
    echo "=== Starting DEV container (live code mount, port 8502) ==="
    mkdir -p output_dev
    docker compose -f docker-compose.dev.yml up -d
    echo "✓ Dev live at https://dev.andom2.bioinfo-wuerz.eu"
    echo "  Edit app.py / src/ directly — changes visible immediately (Streamlit auto-reload)"
    docker compose -f docker-compose.dev.yml logs --tail=5
    ;;

  dev-build)
    echo "=== Rebuilding DEV image (use after requirements.txt change) ==="
    docker compose -f docker-compose.dev.yml build --no-cache
    docker compose -f docker-compose.dev.yml up -d
    echo "✓ Dev rebuilt and running"
    ;;

  dev-stop)
    echo "=== Stopping DEV container ==="
    docker compose -f docker-compose.dev.yml down
    echo "✓ Dev stopped (prod still running)"
    ;;

  status)
    echo "=== Container status ==="
    docker ps --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}" | grep -E "andom2|NAMES"
    echo ""
    echo "Prod: https://andom2.bioinfo-wuerz.eu   (port 8501)"
    echo "Dev:  https://dev.andom2.bioinfo-wuerz.eu (port 8502)"
    ;;

  logs)
    TARGET="${2:-prod}"
    if [ "$TARGET" = "prod" ]; then
      docker compose -f docker-compose.yml logs --tail=30 -f
    else
      docker compose -f docker-compose.dev.yml logs --tail=30 -f
    fi
    ;;

  *)
    echo "Usage: ./deploy.sh [prod|dev|dev-build|dev-stop|status|logs]"
    exit 1
    ;;
esac
