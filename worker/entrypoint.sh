#!/bin/bash
# =============================================================================
# Worker Entrypoint
#
# stanassay is normally baked into the image (compiled at docker build time).
# If you mount a newer /stanassay volume, this will reinstall from source
# on startup — useful for dev iteration without full image rebuild.
# =============================================================================

set -e

echo "============================================"
echo "Immunoplex Batch Worker starting..."
echo "  Redis: ${REDIS_HOST}:${REDIS_PORT}"
echo "  DB:    ${DB_HOST}:${DB_PORT}/${DB_NAME}"
echo "  Output: ${OUTPUT_DIR}"
echo "============================================"

# If /stanassay is mounted (dev override), reinstall from that source.
# This lets you test stanassay changes without rebuilding the whole image.
if [ -d "/stanassay" ] && [ -f "/stanassay/DESCRIPTION" ]; then
  MOUNTED_VER=$(grep '^Version:' /stanassay/DESCRIPTION | awk '{print $2}')
  INSTALLED_VER=$(R --quiet -e 'cat(as.character(packageVersion("stanassay")))' 2>/dev/null || echo "none")
  echo "stanassay mounted: v${MOUNTED_VER} (installed: v${INSTALLED_VER})"

  if [ "$MOUNTED_VER" != "$INSTALLED_VER" ]; then
    echo "Version mismatch — reinstalling stanassay from /stanassay ..."
    echo "  (This compiles Stan models and takes ~5-10 min)"
    R -e "devtools::install('/stanassay', upgrade='never', quiet=FALSE)"
  else
    echo "stanassay versions match — using baked-in build."
  fi
fi

# Verify stanassay is available
if ! R --quiet -e "library(stanassay)" &>/dev/null; then
  echo "ERROR: stanassay is not installed!"
  echo "  Either:"
  echo "    1. Rebuild image with tarball: cp stanassay_*.tar.gz worker/ && docker compose build worker"
  echo "    2. Mount stanassay source: add '- ../stanassay:/stanassay:ro' to docker-compose volumes"
  exit 1
fi

echo "stanassay version: $(R --quiet -e 'cat(as.character(packageVersion("stanassay")))' 2>/dev/null)"
echo "R version: $(R --version | head -1)"
echo "Cores available: $(nproc)"
echo "============================================"

# Create output directory
mkdir -p "${OUTPUT_DIR}"


# Start the Python supervisor
exec python3 /app/supervisor.py
