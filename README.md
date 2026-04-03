# Immunoplex Batch Runner

Generic background job system for running computationally expensive calculations (Bayesian curve fitting, frequentist methods, QC reports, etc.) outside the i-spi Shiny app's R session. Each job specifies a `script_type` that determines which R or Python worker script runs. Jobs are submitted via a REST API, queued in Redis, and processed by worker containers that save results directly to PostgreSQL.

## Architecture

```
                             +------------------+
                             |   i-spi / curl   |
                             |  (job submitter) |
                             +--------+---------+
                                      |
                              POST /jobs (JSON)
                                      |
                             +--------v---------+
                             |   FastAPI (api)   |
                             |    Port 8000      |
                             |  X-API-Key auth   |
                             +--------+---------+
                                      |
                              RPUSH job_id to Redis
                                      |
                             +--------v---------+
                             |     Redis 7       |
                             |  Queue + Status   |
                             +--------+---------+
                                      |
                              BLPOP (blocking pop)
                                      |
                             +--------v---------+
                             |  Python Supervisor |
                             |  (supervisor.py)   |
                             +--------+---------+
                                      |
                              SCRIPT_REGISTRY[script_type] → subprocess
                                      |
                             +--------v---------+
                             |  Worker Script     |
                             |  bayesian → worker_batch.R (stanassay)
                             |  freq     → worker_freq.R  (future)
                             |  custom   → worker_*.py    (any)
                             +--------+---------+
                                      |
                              Upserts per-combo
                                      |
                             +--------v---------+
                             |   PostgreSQL       |
                             |  madi_results.*    |
                             +-------------------+
```

### How It Works

1. **Submit** -- POST a job to the API with `script_type`, study, experiment, scope, and script-specific `params`.
2. **Queue** -- API writes job metadata (including `script_type` and `params` JSON) to a Redis hash and pushes the job ID onto a Redis list.
3. **Pick up** -- Worker's Python supervisor BLPOPs the queue, reads job metadata, looks up `script_type` in `SCRIPT_REGISTRY` to find the right interpreter + script.
4. **Dispatch** -- Supervisor spawns the script with common CLI args (`--study`, `--scope`, etc.) plus script-specific args from `params` dict (each key becomes `--key value`).
5. **Execute** -- The script does its work (e.g. Bayesian fitting, frequentist curves, QC) and saves results to PostgreSQL after each combo.
6. **Monitor** -- Supervisor reads the progress file every 5 seconds and pushes percentage, ETA, and current antigen info to Redis.
7. **Poll** -- Client polls `GET /jobs/{id}` to see real-time progress, or waits for `status: "completed"`.

### Key Design Decisions

- **Script registry**: Adding a new calculation is just dropping an R/Python file and adding one line to `SCRIPT_REGISTRY` in `supervisor.py`. No API changes needed.
- **Generic params**: Script-specific arguments are passed via a `params` dict, converted to `--key value` CLI args. Each script takes what it needs, ignores the rest.
- **Sequential processing**: Combos are processed one at a time (not parallel) so that progress updates and DB saves happen after each combo.
- **DB-first output**: No RData files. All results go directly to PostgreSQL for immediate query by i-spi.
- **Idempotent upserts**: `ON CONFLICT DO UPDATE` ensures re-running a job safely overwrites previous results.
- **Workspace isolation**: `project_id` is part of every natural key, so the same study name in different workspaces stays separate.
- **Coexistence**: Bayesian results (`bayes_*` tables) and frequentist results (`best_*` tables) are completely independent -- both can exist for the same experiment.

---

## Quick Start (Local Development)

### Prerequisites

- Docker Desktop (macOS/Windows) or Docker Engine (Linux)
- PostgreSQL running on host with `local_madi_ispi` database
- `stanassay` R package tarball in `worker/` directory

### 1. Build the stanassay tarball (if not already done)

```bash
cd ../stanassay
R CMD build . --no-manual --no-vignettes
cp stanassay_*.tar.gz ../immunoplex_batch_calculator/worker/
```

### 2. Run the database migration (one-time)

```bash
psql -d local_madi_ispi -f migrations/001_create_bayes_tables.sql
```

### 3. Start the stack

```bash
cd immunoplex_batch_calculator
docker compose up --build
```

This starts 3 containers:
- **redis** -- Redis 7 on host port 6380 (password: `devpassword`)
- **api** -- FastAPI on host port 8000
- **worker** -- R + Python supervisor

First build takes ~10 minutes (Stan C++ compilation). Subsequent builds are fast (cached layers).

### 4. Verify

```bash
# Health check (no auth)
curl http://localhost:8000/health

# Swagger UI
open http://localhost:8000/docs
```

---

## API Reference

All endpoints except `/health` require the `X-API-Key` header.

**Local dev key**: `dev-key-immunoplex`

### POST /jobs -- Submit a Job

```bash
curl -X POST http://localhost:8000/jobs \
  -H "Content-Type: application/json" \
  -H "X-API-Key: dev-key-immunoplex" \
  -d '{
    "project_id": 17,
    "study": "MADI_P3_GAPS",
    "scope": "study"
  }'
```

#### Request Body

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `project_id` | int | **yes** | -- | Workspace/project ID. Critical for isolation. |
| `study` | string | **yes** | -- | Study accession (e.g. `"MADI_P3_GAPS"`) |
| `experiment` | string | no | null | Experiment accession. Required if scope is `"experiment"` or `"antigen"`. |
| `antigen` | string | no | null | Antigen name. Required if scope is `"antigen"`. |
| `source` | string | no | null | Standard source filter (e.g. `"NIBSC06_140"`). If omitted, all sources. |
| `scope` | string | no | `"study"` | `"study"` (all experiments), `"experiment"` (single), or `"antigen"` (single antigen). |
| `script_type` | string | no | `"bayesian"` | Which worker script to run. See `SCRIPT_REGISTRY` in `supervisor.py`. |
| `params` | dict | no | `{}` | Script-specific params. Each key becomes a `--key value` CLI arg. |
| `cdan_cv_threshold` | float | no | `20.0` | Shorthand for bayesian CDAN CV%. Auto-merged into `params` as `cdan_cv`. |

#### Scope Examples

```bash
# Full study -- fits ALL experiments x antigens x plates
curl -X POST http://localhost:8000/jobs \
  -H "Content-Type: application/json" \
  -H "X-API-Key: dev-key-immunoplex" \
  -d '{"project_id":17, "study":"MADI_P3_GAPS", "scope":"study"}'

# Single experiment
curl -X POST http://localhost:8000/jobs \
  -H "Content-Type: application/json" \
  -H "X-API-Key: dev-key-immunoplex" \
  -d '{"project_id":17, "study":"MADI_P3_GAPS", "experiment":"ADCD", "scope":"experiment"}'

# Single antigen (fastest for testing)
curl -X POST http://localhost:8000/jobs \
  -H "Content-Type: application/json" \
  -H "X-API-Key: dev-key-immunoplex" \
  -d '{"project_id":17, "study":"MADI_P3_GAPS", "experiment":"ADCD", "antigen":"tt", "scope":"antigen"}'

# Custom CDAN threshold (25% CV)
curl -X POST http://localhost:8000/jobs \
  -H "Content-Type: application/json" \
  -H "X-API-Key: dev-key-immunoplex" \
  -d '{"project_id":17, "study":"MADI_P3_GAPS", "experiment":"ADCD", "scope":"experiment", "cdan_cv_threshold": 25.0}'
```

### GET /jobs/{job_id} -- Poll Job Status

```bash
curl -H "X-API-Key: dev-key-immunoplex" \
  http://localhost:8000/jobs/{job_id}
```

#### Response

```json
{
  "job_id": "cacebb23-af58-4089-98fe-bb503bd05d08",
  "status": "running",
  "script_type": "bayesian",
  "project_id": 17,
  "study": "MADI_P3_GAPS",
  "experiment": "ADCD",
  "source": null,
  "scope": "antigen",
  "cdan_cv_threshold": 20.0,
  "params": {"cdan_cv": 20.0},
  "created_at": "2026-04-02T23:02:55.881117+00:00",
  "started_at": "2026-04-02T23:02:55.884920+00:00",
  "completed_at": null,
  "progress": "1/2",
  "total_combos": 2,
  "completed_combos": 1,
  "percentage": 50.0,
  "experiment_progress": "0/1",
  "experiments_done": 0,
  "experiments_total": 1,
  "elapsed_minutes": 2.5,
  "eta_minutes": 2.5,
  "eta_display": "~2 min 30 sec",
  "speed_seconds_per_combo": 150.3,
  "current_experiment": "ADCD",
  "current_antigens": "tt",
  "output_path": null,
  "error": null
}
```

**Status values**: `queued` -> `running` -> `completed` | `failed` | `cancelled`

**Progress fields for UI**:

| Field | Description |
|-------|-------------|
| `percentage` | 0.0 - 100.0, suitable for a progress bar |
| `progress` | Human-readable combo progress like `"3/10"` |
| `experiment_progress` | Experiment-level like `"1/4"` |
| `eta_display` | Human-readable time remaining like `"~3 min 20 sec"` |
| `eta_minutes` | Numeric minutes remaining (for calculations) |
| `speed_seconds_per_combo` | Average seconds per combo (for display) |
| `elapsed_minutes` | Wall-clock time since job started |
| `current_experiment` | Which experiment is being processed now |
| `current_antigens` | Which antigen is being fitted now |

### GET /jobs -- List All Jobs

```bash
# All jobs (newest first)
curl -H "X-API-Key: dev-key-immunoplex" \
  http://localhost:8000/jobs

# Filter by study
curl -H "X-API-Key: dev-key-immunoplex" \
  "http://localhost:8000/jobs?study=MADI_P3_GAPS"

# Filter by status
curl -H "X-API-Key: dev-key-immunoplex" \
  "http://localhost:8000/jobs?status=running"

# Filter by script type
curl -H "X-API-Key: dev-key-immunoplex" \
  "http://localhost:8000/jobs?script_type=bayesian"

# Pagination
curl -H "X-API-Key: dev-key-immunoplex" \
  "http://localhost:8000/jobs?limit=10&offset=20"
```

### DELETE /jobs/{job_id} -- Cancel a Job

```bash
curl -X DELETE -H "X-API-Key: dev-key-immunoplex" \
  http://localhost:8000/jobs/{job_id}
```

If the job is **queued**, the worker will skip it. If **running**, the worker terminates the R process after the current combo.

### GET /health -- Health Check

```bash
curl http://localhost:8000/health
# {"status":"ok","redis":"connected"}
```

No authentication required.

---

## Database Schema

Results are stored in the `madi_results` schema across 3 tables. See `migrations/001_create_bayes_tables.sql` for full DDL.

### `bayes_curves` -- Curve parameters + diagnostics

One row per **plate x antigen x source** combo. Contains:
- Posterior median curve params: `a`, `b`, `c`, `d`, `g`
- Model selection: `curve_family`, `plate_best_family`, `global_best_family`
- CDAN LOQ at 20%: `lloq`, `uloq` (with y-values)
- CDAN LOQ at 15%: `lloq_15`, `uloq_15`
- CDAN LOQ at custom threshold: `lloq_custom`, `uloq_custom`, `cdan_cv_threshold`
- Second derivative limits: `lo2d`, `uo2d`
- Detection limits: `lod`, `lrdl`, `uod`, `urdl`
- Inflection point with 95% CIs: `inflect_x`, `inflect_y`, `inflect_x_lower`, `inflect_x_upper`
- Per-plate ELPD for each family, global stacking weights

**Natural key**: `(project_id, study_accession, experiment_accession, plateid, plate, nominal_sample_dilution, source, wavelength, antigen, feature)`

### `bayes_samples` -- Sample concentrations

One row per **sample well x dilution**. Contains:
- `raw_predicted_concentration` (posterior median)
- `conc_lower`, `conc_upper` (95% credible interval)
- `se_concentration`, `pcov`
- `gate_class` (in range, below LOQ, etc.)
- FK: `bayes_curves_id`

### `bayes_ensemble` -- Per-family model comparison

One row per **plate x antigen x family** (3 rows per plate: 4pl, 5pl, gompertz). Contains:
- `plate_elpd`, `is_plate_best`
- `global_stacking_weight`, `is_global_best`
- Parameter posteriors with 95% CIs: `a`, `a_lower`, `a_upper`, etc.
- FK: `bayes_curves_id`

### Querying Results

```sql
-- All curves for a study
SELECT experiment_accession, antigen, plateid, curve_family, lloq, uloq
FROM madi_results.bayes_curves
WHERE project_id = 17 AND study_accession = 'MADI_P3_GAPS';

-- Sample concentrations for a specific antigen
SELECT s.sampleid, s.dilution, s.mfi, s.raw_predicted_concentration,
       s.conc_lower, s.conc_upper, s.gate_class
FROM madi_results.bayes_samples s
JOIN madi_results.bayes_curves c ON c.bayes_curves_id = s.bayes_curves_id
WHERE c.study_accession = 'MADI_P3_GAPS'
  AND c.experiment_accession = 'ADCD'
  AND c.antigen = 'tt';

-- Model comparison across families
SELECT plateid, family, plate_elpd, is_plate_best, global_stacking_weight
FROM madi_results.bayes_ensemble
WHERE study_accession = 'MADI_P3_GAPS' AND experiment_accession = 'ADCD'
ORDER BY plateid, family;
```

### Reconstructing Plots from DB

The `bayes_curves` table stores everything needed to reconstruct `stanassay::plot_bayesian_plate()`:
- **Fitted curve**: evaluate `f(x) = d + (a - d) / (1 + (x/c)^b)^g` using stored params
- **Standards overlay**: query `xmap_standard` for the same plate (already in DB)
- **Sample overlay**: query `bayes_samples` for MFI vs concentration points
- **CI ribbons**: use `bayes_ensemble` per-family param CIs
- **LOQ lines**: vertical lines from `lloq`, `uloq`, `lloq_custom`, `uloq_custom`
- **LOD/inflection**: from stored values in `bayes_curves`

---

## Project Structure

```
immunoplex_batch_calculator/
  api/
    app.py              # FastAPI application
    Dockerfile          # Python 3.11-slim
    requirements.txt    # fastapi, uvicorn, redis
  worker/
    worker_batch.R      # R worker script (stanassay fitting + DB save)
    supervisor.py       # Python supervisor (BLPOP + subprocess + progress)
    entrypoint.sh       # Container entrypoint (stanassay version check)
    Dockerfile          # rocker/tidyverse + rstan + stanassay + Python
    requirements.txt    # redis
    stanassay_*.tar.gz  # Private package tarball (not in git)
  migrations/
    001_create_bayes_tables.sql   # One-time DB migration
  k8s/
    immunoplex-batch-cal-api.yml     # K8s Deployment + Service
    immunoplex-batch-cal-worker.yml  # K8s Deployment
  docker-compose.yml    # Local dev stack
  .env.example          # Environment variable reference
```

---

## Adding a New Calculation Script

The worker is a **generic batch runner**. The supervisor dispatches jobs to scripts based on `script_type`, using a registry in `supervisor.py`:

```python
# worker/supervisor.py — SCRIPT_REGISTRY
SCRIPT_REGISTRY = {
    "bayesian":    ("Rscript", SCRIPTS_DIR / "worker_batch.R"),
    # "frequentist": ("Rscript", SCRIPTS_DIR / "worker_freq.R"),
    # "qc_report":   ("python3", SCRIPTS_DIR / "worker_qc.py"),
}
```

To add a new script (e.g. frequentist), you only need to:
1. Write the script (R or Python)
2. Add one line to `SCRIPT_REGISTRY`
3. Add `COPY` + package installs to `Dockerfile`
4. Submit jobs with `"script_type": "frequentist"`

No API changes, no queue changes, no supervisor logic changes.

### Step 1: Write the R Script

Create `worker/worker_freq.R` following this contract:

```r
# =============================================================================
# worker_freq.R -- Frequentist Batch Worker
#
# CLI args: --study, --experiment, --antigen, --source, --scope,
#           --job_id, --project_id, --output_dir, --progress_dir
# DB creds: from env vars DB_HOST, DB_PORT, DB_NAME, DB_USER, DB_PASSWORD
# Progress: write JSON to {progress_dir}/progress_{job_id}.json
# =============================================================================

suppressPackageStartupMessages({
  library(RPostgres)
  library(DBI)
  library(dplyr)
  library(jsonlite)
  # library(your_package)
})

# ---- CLI Arguments ----
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list(
    study = "", experiment = "", antigen = "", source = "",
    scope = "study", job_id = "local", project_id = "0",
    output_dir = "batch_output", progress_dir = "/tmp"
    # Add any script-specific params here
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    param_name <- sub("^--", "", key)
    if (param_name %in% names(params)) {
      params[[param_name]] <- args[i + 1L]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  params
}
PARAMS <- parse_args()

# ---- DB Connection ----
open_conn <- function() {
  DBI::dbConnect(RPostgres::Postgres(),
    dbname   = Sys.getenv("DB_NAME", "local_madi_ispi"),
    host     = Sys.getenv("DB_HOST", "localhost"),
    port     = as.integer(Sys.getenv("DB_PORT", "5432")),
    user     = Sys.getenv("DB_USER", "hardik"),
    password = Sys.getenv("DB_PASSWORD", "")
  )
}

# ---- Progress Reporting ----
# The supervisor reads this file and pushes updates to Redis.
# Required fields: total_combos, completed_combos
# Optional: current_experiment, current_antigens, experiments_done, experiments_total
write_progress <- function(total, completed, experiment, antigen,
                           exp_done = 0, exp_total = 0) {
  jsonlite::write_json(
    list(
      total_combos = total,
      completed_combos = completed,
      current_experiment = experiment,
      current_antigens = antigen,
      experiments_done = exp_done,
      experiments_total = exp_total
    ),
    file.path(PARAMS$progress_dir, paste0("progress_", PARAMS$job_id, ".json")),
    auto_unbox = TRUE
  )
}

# ---- Main Loop ----
conn <- open_conn()

# Discover combos (your logic here)
combos <- discover_combos(conn, PARAMS$study, PARAMS$project_id, ...)

total <- nrow(combos)
completed <- 0L
write_progress(total, completed, "", "")

for (i in seq_len(total)) {
  combo <- combos[i, ]

  # Your fitting logic
  result <- fit_frequentist(conn, combo)

  # Save to DB (your upsert logic)
  save_to_db(conn, result, PARAMS$job_id)

  completed <- completed + 1L
  write_progress(total, completed, combo$experiment, combo$antigen)
}

DBI::dbDisconnect(conn)
cat(sprintf("\n==== BATCH COMPLETE -- %d/%d combos ====\n", completed, total))
```

**Key contract points**:
1. Accept `--job_id` and `--progress_dir` CLI args
2. Write a `progress_{job_id}.json` file with at least `total_combos` and `completed_combos`
3. Read DB creds from environment variables
4. Print to stdout (supervisor captures it for docker logs)
5. Exit 0 on success, non-zero on failure

### Step 2: Register in SCRIPT_REGISTRY

Add one line to `worker/supervisor.py`:

```python
SCRIPT_REGISTRY = {
    "bayesian":     ("Rscript", SCRIPTS_DIR / "worker_batch.R"),
    "frequentist":  ("Rscript", SCRIPTS_DIR / "worker_freq.R"),  # ← add this
}
```

That's it. The supervisor already handles:
- Looking up the interpreter + script path
- Passing common args (`--study`, `--scope`, `--job_id`, etc.)
- Passing script-specific args from `params` dict as `--key value`
- Progress monitoring, cancellation, error capture

### Step 3: Add Dependencies to Dockerfile

If your new script needs additional R packages, add them to `worker/Dockerfile`:

```dockerfile
# Additional R packages for frequentist worker
RUN R -e "install.packages(c('nplr', 'minpack.lm'), repos='https://cloud.r-project.org/')"
```

And copy the script:

```dockerfile
COPY worker_freq.R .
```

### Step 4: Create DB Tables (if needed)

Add a new migration in `migrations/`:

```
migrations/002_create_freq_tables.sql
```

Run it once:

```bash
psql -d local_madi_ispi -f migrations/002_create_freq_tables.sql
```

### Step 5: Submit Jobs with the New Type

```bash
# Frequentist job with script-specific params
curl -X POST http://localhost:8000/jobs \
  -H "Content-Type: application/json" \
  -H "X-API-Key: dev-key-immunoplex" \
  -d '{
    "project_id": 17,
    "study": "MADI_P3_GAPS",
    "experiment": "ADCD",
    "scope": "experiment",
    "script_type": "frequentist",
    "params": {"method": "nplr", "weights": "1/y2"}
  }'
```

The `params` dict values are passed to the R script as `--method nplr --weights 1/y2`.

### Adding a Python Script Instead of R

The registry supports any interpreter. Python scripts work the same way:

```python
SCRIPT_REGISTRY = {
    "bayesian":     ("Rscript", SCRIPTS_DIR / "worker_batch.R"),
    "frequentist":  ("Rscript", SCRIPTS_DIR / "worker_freq.R"),
    "qc_report":    ("python3", SCRIPTS_DIR / "worker_qc.py"),
}
```

The only requirement is that the script follows the progress file contract (write `progress_{job_id}.json` with `total_combos` and `completed_combos`).

---

## Docker Operations

```bash
# Start everything
docker compose up --build

# Rebuild only the worker (after changing R code)
docker compose build worker && docker compose up -d worker

# Rebuild without cache (when Docker gets confused)
docker compose build --no-cache worker

# View worker logs (R output prefixed with [R])
docker compose logs -f worker

# Stop everything (clears Redis queue)
docker compose down

# Stop without losing Redis data
docker compose stop
```

### Updating stanassay

When the stanassay R package changes:

```bash
cd ../stanassay
R CMD build . --no-manual --no-vignettes
cp stanassay_*.tar.gz ../immunoplex_batch_calculator/worker/
cd ../immunoplex_batch_calculator
docker compose build worker   # ~5-10 min (Stan recompilation)
docker compose up -d worker
```

**Dev shortcut**: Mount stanassay source to skip image rebuild (slower startup but faster iteration):

```yaml
# In docker-compose.yml, uncomment:
volumes:
  - ../stanassay:/stanassay:ro
```

The entrypoint will detect the mount and reinstall from source if the version differs.

---

## Kubernetes Deployment

Manifests are in `k8s/`. The worker needs significant resources for Stan MCMC.

```bash
# Apply (after building and pushing images to your registry)
kubectl apply -f k8s/immunoplex-batch-cal-api.yml
kubectl apply -f k8s/immunoplex-batch-cal-worker.yml
```

Key differences from local dev:
- Redis and PostgreSQL connection details come from Kubernetes secrets (`madi` sealed secret)
- Worker connects to the prod/preprod PostgreSQL instance
- API is exposed as a ClusterIP service (accessed via ingress or port-forward)
- Worker resource limits: 4-8 CPU, 8-16Gi RAM

---

## Environment Variables

| Variable | Default | Used By | Description |
|----------|---------|---------|-------------|
| `REDIS_HOST` | `redis` | api, worker | Redis hostname |
| `REDIS_PORT` | `6379` | api, worker | Redis port |
| `REDIS_AUTH` | (empty) | api, worker | Redis password |
| `REDIS_DB` | `0` | api, worker | Redis database number |
| `API_KEY` | `dev-key-immunoplex` | api | API authentication key |
| `DB_NAME` | `local_madi_ispi` | worker | PostgreSQL database |
| `DB_HOST` | `localhost` | worker | PostgreSQL host |
| `DB_PORT` | `5432` | worker | PostgreSQL port |
| `DB_USER` | `hardik` | worker | PostgreSQL user |
| `DB_PASSWORD` | (empty) | worker | PostgreSQL password |
| `DB_SSLMODE` | `disable` | worker | PostgreSQL SSL mode |
| `OUTPUT_DIR` | `/data/bayes` | worker | Output directory (legacy) |
| `PROGRESS_DIR` | `/tmp` | worker | Progress file directory |

---

## Troubleshooting

**Job stuck in "queued"**: Worker might not be running. Check `docker compose logs worker`.

**Job fails immediately**: Check the `error` field in job status. Common causes:
- PostgreSQL not reachable from worker container (check `DB_HOST`)
- Study/experiment not found in database
- stanassay not installed (check entrypoint logs)

**Stan warnings in logs**: Divergent transitions and max treedepth warnings are normal for some plates. The ensemble model selection handles this by downweighting poorly-fitting families.

**Docker build uses old code**: Use `docker compose build --no-cache worker` to force a clean rebuild.

**Stale jobs in queue after rebuild**: `docker compose down` (not just `stop`) clears Redis.

**DB connection from worker**: The worker container uses `host.docker.internal` to reach the host's PostgreSQL. On Linux, you may need `--add-host=host.docker.internal:host-gateway` in the compose file.
