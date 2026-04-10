# =============================================================================
# Worker Dockerfile — R + stanassay + Python supervisor
#
# Base: rocker/tidyverse (includes R, tidyverse, devtools)
# Adds: rstan, stanassay (from GitHub), Python 3, redis-py
# Entrypoint: Python supervisor that BLPOP on Redis and spawns Rscript
#
# stanassay is installed from GitHub (immunoplex/stanassay).
# The Stan C++ compilation happens ONCE at build time and is cached in the
# Docker layer — container startup is instant.
#
# Build:
#   docker compose build worker        # uses GitHub install
#   docker compose up --build           # build + run
#
# Local dev: install directly from a local clone if needed:
#   R -e "devtools::install_local('../stanassay')"
# =============================================================================

FROM rocker/tidyverse:latest

LABEL org.opencontainers.image.source=https://github.com/immunoplex/immunoplex-batch-cal-worker
LABEL org.opencontainers.image.description="Immunoplex Batch Runner — Worker"

# ── System dependencies ──────────────────────────────────────────────────────

RUN apt-get update && apt-get install -y --no-install-recommends \
    libpq-dev \
    python3 \
    python3-pip \
    python3-venv \
    zip \
    && rm -rf /var/lib/apt/lists/*

# ── Python dependencies (in a venv to avoid externally-managed-environment) ──

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
COPY requirements.txt /tmp/py-requirements.txt
RUN pip install --no-cache-dir -r /tmp/py-requirements.txt

# ── R dependencies for stanassay ─────────────────────────────────────────────

# Install rstan and Stan toolchain first (takes a while but caches well)
RUN R -e "install.packages(c('rstan', 'rstantools', 'StanHeaders', 'loo', 'R6', 'BH', 'Rcpp', 'RcppEigen', 'RcppParallel'), repos='https://cloud.r-project.org/')"

# Additional R deps for the worker script
RUN R -e "install.packages(c('RPostgres', 'DBI', 'dplyr', 'ggplot2', 'plotly', 'jsonlite'), repos='https://cloud.r-project.org/')"

# ── Install stanassay from GitHub ────────────────────────────────────────────
# stanassay is a public package at immunoplex/stanassay.
# Stan C++ → .so compilation happens here at build time (~5-10 min, cached).
RUN R -e "devtools::install_github('immunoplex/stanassay', upgrade='never')"

# ── Copy worker scripts ─────────────────────────────────────────────────────

WORKDIR /app
COPY supervisor.py .
COPY worker_batch.R .
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh

# ── Environment defaults ─────────────────────────────────────────────────────

# Runtime defaults — credentials are set via K8s secrets / docker-compose env
ENV REDIS_HOST=redis \
    REDIS_PORT=6379 \
    REDIS_DB=0 \
    DB_NAME=local_madi_ispi \
    DB_HOST=localhost \
    DB_PORT=5432 \
    DB_USER=postgres \
    DB_SSLMODE=disable \
    OUTPUT_DIR=/data/bayes \
    PROGRESS_DIR=/tmp

# ── Entrypoint ───────────────────────────────────────────────────────────────

CMD ["/app/entrypoint.sh"]
