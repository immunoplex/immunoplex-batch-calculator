# =============================================================================
# worker_batch.R — Bayesian Standard Curve Batch Fitting Worker
#
# Adapted from i-spi/research/hpc_bayes_batch.R and i-spi/src/std_curver_ui.R
# for use in the immunoplex_batch_calculator worker container.
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  OVERVIEW                                                                   │
# │                                                                             │
# │  This script fits Bayesian hierarchical standard curves (4PL, 5PL,         │
# │  Gompertz ensemble) to immunoassay data using the stanassay R package.     │
# │  It is spawned as a subprocess by supervisor.py when a "bayesian" job      │
# │  is dequeued from Redis.                                                    │
# │                                                                             │
# │  For each antigen × source combo:                                           │
# │    1. Fetches standards, samples, and blanks from PostgreSQL                │
# │    2. Converts raw dilutions to concentrations using the correct scale      │
# │    3. Applies prozone correction to standards                               │
# │    4. Fits a 3-family ensemble (4PL, 5PL, Gompertz) via Stan MCMC         │
# │    5. Computes diagnostics (LOQ, LOD, LRDL, URDL, inflection, etc.)       │
# │    6. Back-calculates sample concentrations with credible intervals         │
# │    7. Saves everything to PostgreSQL (5 tables) after each combo           │
# │    8. Writes a progress JSON file for real-time monitoring                  │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  CONCENTRATION SCALING — CRITICAL                                           │
# │                                                                             │
# │  Immunoassay standards are measured at serial dilutions (e.g. 10, 20, 40,  │
# │  ..., 2560). The MFI (signal) is recorded at each dilution. To convert     │
# │  dilution factors to concentration (the x-axis for curve fitting), we use: │
# │                                                                             │
# │     concentration = base_num / dilution_factor                              │
# │                                                                             │
# │  Where base_num is the UNDILUTED concentration of the standard. This is    │
# │  stored in madi_results.xmap_antigen_family.standard_curve_concentration.  │
# │                                                                             │
# │  Example for antigen "tt" in ADCD:                                          │
# │    standard_curve_concentration = 10000 (IU/mL or ng/mL)                   │
# │    dilution=10 → conc = 10000/10 = 1000                                    │
# │    dilution=2560 → conc = 10000/2560 = 3.9                                │
# │    → x-axis range: log10(3.9) ≈ 0.6 to log10(1000) ≈ 3.0                 │
# │                                                                             │
# │  FALLBACK: If standard_curve_concentration is not found in the antigen     │
# │  family table, falls back to nominal_sample_dilution (from xmap_header).   │
# │  This is typically wrong (10 instead of 10000) and will produce curves     │
# │  on the wrong scale. Always ensure the antigen family table is populated.  │
# │                                                                             │
# │  This matches i-spi's behavior at std_curver_ui.R line 2620:              │
# │    base_num = sc_concentration_cache() %||% coalesce(nom, 100000)          │
# │  where sc_concentration_cache() is the standard_curve_concentration.       │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  DATA PIPELINE                                                              │
# │                                                                             │
# │  1. COMBO DISCOVERY (discover_combos)                                       │
# │     Query: JOIN xmap_standard × xmap_header                                │
# │     Filters: project_id, study, experiment, antigen, source                │
# │     Result: unique (experiment, antigen, source, feature, wavelength) tuples│
# │     Each tuple is one "combo" = one independent curve fit.                  │
# │                                                                             │
# │  2. FETCH STANDARDS (fetch_standards)                                       │
# │     Query: xmap_standard JOIN xmap_header                                  │
# │     Filtered by: study, experiment, antigen, wavelength, source            │
# │     Returns: plateid, dilution_factor, mfi, nominal_sample_dilution, etc.  │
# │     IMPORTANT: Source filter means each source (e.g. NIBSC06_140 vs SD)    │
# │     is fitted INDEPENDENTLY. Same plate can have 2 fits (one per source).  │
# │                                                                             │
# │  3. CONCENTRATION CONVERSION                                                │
# │     base_num = standard_curve_concentration from xmap_antigen_family       │
# │     concentration = base_num / dilution_factor                              │
# │     Then: group_by(plateid, concentration) → median(mfi) per point         │
# │     This aggregates duplicate standard wells (same dilution) per plate.    │
# │                                                                             │
# │  4. PROZONE CORRECTION (correct_prozone)                                    │
# │     High-dose hook effect: at very high concentrations, some assays show   │
# │     DECREASING MFI. This confuses curve fitting (multiple x for same y).   │
# │     Fix: for standards beyond the MFI peak, adjust MFI upward using a      │
# │     dampened correction formula. Parameters: prop_diff=0.1, dil_scale=2.   │
# │     Applied per-plate BEFORE fitting.                                       │
# │     Inlined from i-spi's prozone correction (not stanassay's).             │
# │                                                                             │
# │  5. STANASSAY ENSEMBLE FIT (StanAssay$fit_ensemble)                         │
# │     Creates StanAssay R6 object with:                                       │
# │       std_data = prozone-corrected standards (all plates pooled)           │
# │       concentration_col = "concentration"                                   │
# │       response_col = "mfi"                                                  │
# │       plate_col = "plateid" (hierarchical: plates share priors)            │
# │       blank_data = buffer well MFI (for LOD computation)                   │
# │       apply_prozone = FALSE (already applied above)                        │
# │     Fits 3 families: 4PL, 5PL, Gompertz                                   │
# │     Stan MCMC: 4 chains × 1000 iterations = 2000 posterior draws/chain    │
# │     Model selection: LOO-ELPD stacking weights across families             │
# │     Result: per-plate best family + global best family + stacking weights  │
# │                                                                             │
# │  6. PARAMETER EXTRACTION (extract_plate_params_v2)                          │
# │     For each plate, extracts posterior MEDIAN of (a, b, c, d, g) from the  │
# │     plate's best-family fit. These are the "point estimate" curve params.  │
# │     Note: c is stored in linear scale (exp(log_c)) — Stan fits log_c.     │
# │                                                                             │
# │  7. DIAGNOSTICS (per plate)                                                 │
# │     All computed by stanassay's methods:                                    │
# │     - CDAN precision profile: CV% vs concentration via delta method         │
# │       → LLOQ/ULOQ at 15%, 20%, and custom CV thresholds                   │
# │     - LOD: Bayesian limit of detection (3σ blank criterion)                │
# │     - LRDL/URDL: posterior predictive reliable detection limits            │
# │     - LO2D/UO2D: second-derivative shoulder points (exact closed-form)    │
# │     - Inflection: max |dy/dx| point with 95% credible interval            │
# │                                                                             │
# │  8. SAMPLE BACK-CALCULATION (compute_se_pcov_v2)                           │
# │     For each sample MFI value:                                              │
# │     a) Use stanassay's predict_concentration_bayesian() for point est + CI │
# │     b) Extract full posterior via curve inverse: x = f⁻¹(MFI)             │
# │        → raw_predicted_concentration = median of posterior draws           │
# │        → se_concentration = SD of posterior draws                           │
# │        → pcov = SD/median (coefficient of variation)                       │
# │        → conc_lower/conc_upper = 2.5%/97.5% credible interval             │
# │     c) Gate classification: Below LLOQ / Within Range / Above ULOQ        │
# │     Concentrations are in the same units as the standard curve x-axis      │
# │     (i.e. using standard_curve_concentration scaling).                      │
# │                                                                             │
# │  9. PLOT DATA (curve_grid + cdan_grid)                                     │
# │     Pre-computed grids for reconstructing the interactive plot from DB:     │
# │     - curve_grid: 100 log-spaced x-points with posterior median MFI        │
# │       and 95% CI bounds. Uses 200 correlated MCMC draws (NOT marginal CIs)│
# │       via stanassay::build_plot_data(). Stored in bayes_curve_grid.        │
# │     - cdan_grid: 200-point CDAN precision profile (smoothed CV% vs conc). │
# │       Uses the full noise model (theta_base, theta_prop, gamma) +          │
# │       analytical delta method + LOESS smoothing.                            │
# │       Stored in bayes_cdan_grid.                                            │
# │     These grids make it possible to reconstruct plot_bayesian_plate()      │
# │     output from DB alone, without needing the StanAssay fit object.        │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  DATABASE OUTPUT (5 tables in madi_results schema)                          │
# │                                                                             │
# │  bayes_curves      1 row per plate × source                                │
# │    Curve params (a,b,c,d,g), LOQ at 15/20/custom%, LOD, LRDL, URDL,       │
# │    inflection, ELPD, stacking weights. Parent table for FK.                │
# │                                                                             │
# │  bayes_samples     1 row per sample well × dilution                        │
# │    MFI, predicted concentration (median + 95% CI), SE, pCoV, gate_class.   │
# │    FK → bayes_curves.                                                       │
# │                                                                             │
# │  bayes_ensemble    3 rows per plate (one per family: 4pl, 5pl, gompertz)   │
# │    Per-family ELPD, stacking weight, parameter posteriors with 95% CIs.    │
# │    FK → bayes_curves.                                                       │
# │                                                                             │
# │  bayes_curve_grid  100 rows per plate × source                             │
# │    Pre-computed fitted curve + 95% CI ribbon for plot reconstruction.       │
# │    FK → bayes_curves.                                                       │
# │                                                                             │
# │  bayes_cdan_grid   200 rows per plate × source                             │
# │    Pre-computed CDAN precision profile for plot reconstruction.             │
# │    FK → bayes_curves.                                                       │
# │                                                                             │
# │  All writes use idempotent upsert (INSERT ... ON CONFLICT DO UPDATE)       │
# │  or delete+insert (for grid tables). Re-running a job safely overwrites.   │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  CLI ARGUMENTS (passed by supervisor.py)                                    │
# │                                                                             │
# │  --study          Study accession (e.g. MADI_P3_GAPS) — required           │
# │  --experiment     Experiment accession (e.g. ADCD) — for experiment/antigen│
# │  --antigen        Antigen name (e.g. tt) — for antigen scope               │
# │  --source         Standard source filter (e.g. NIBSC06_140)                │
# │  --scope          "study" | "experiment" | "antigen"                        │
# │  --job_id         UUID for tracking in Redis and DB                         │
# │  --project_id     Workspace ID (for multi-tenant isolation)                │
# │  --cdan_cv        Custom CDAN CV% threshold (default 20)                   │
# │  --output_dir     Legacy (not used — all output goes to DB)                │
# │  --progress_dir   Where to write progress_{job_id}.json                    │
# │                                                                             │
# │  Additional --key value pairs from the API's params dict are accepted      │
# │  (the arg parser is permissive) but currently unused by this script.       │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  ENVIRONMENT VARIABLES (DB credentials, set by Docker)                      │
# │                                                                             │
# │  DB_NAME, DB_HOST, DB_PORT, DB_USER, DB_PASSWORD, DB_SSLMODE              │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# ┌─────────────────────────────────────────────────────────────────────────────┐
# │  PROCESSING ORDER                                                           │
# │                                                                             │
# │  Sequential by design (not parallel). For each combo:                       │
# │    fit_one() → save_to_db() → write_progress()                             │
# │  This ensures:                                                              │
# │    - Real-time progress reporting after each combo                          │
# │    - Partial results available in DB if the job is cancelled/fails          │
# │    - Memory doesn't accumulate across combos                                │
# └─────────────────────────────────────────────────────────────────────────────┘
#
# Dependencies: stanassay, RPostgres, DBI, dplyr, parallel, jsonlite
# =============================================================================


suppressPackageStartupMessages({
  library(parallel)
  library(RPostgres)
  library(DBI)
  library(dplyr)
  library(stanassay)
  library(jsonlite)
})


# ─── CLI Arguments ───────────────────────────────────────────────────────────

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults
  params <- list(
    study        = "",
    experiment   = "",
    antigen      = "",
    source       = "",
    scope        = "study",
    job_id       = "local",
    project_id   = "0",
    cdan_cv      = "20",
    output_dir   = "batch_output",
    progress_dir = "/tmp"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (grepl("^--", key)) {
      param_name <- sub("^--", "", key)
      if (i + 1L <= length(args)) {
        params[[param_name]] <- args[i + 1L]
        i <- i + 2L
      } else {
        stop(sprintf("Missing value for %s", key))
      }
    } else {
      warning(sprintf("Ignoring positional argument: %s", key))
      i <- i + 1L
    }
  }

  # Convert numeric params
  params$project_id <- as.integer(params$project_id)
  params$cdan_cv    <- as.numeric(params$cdan_cv)

  if (nchar(params$study) == 0L) stop("--study is required")
  if (params$scope %in% c("experiment", "antigen") && nchar(params$experiment) == 0L) {
    stop("--experiment is required when scope='experiment' or 'antigen'")
  }
  if (params$scope == "antigen" && nchar(params$antigen) == 0L) {
    stop("--antigen is required when scope='antigen'")
  }

  params
}


PARAMS <- parse_args()

message(sprintf("Worker starting: study=%s experiment=%s scope=%s job_id=%s project_id=%d cdan_cv=%.1f",
                PARAMS$study, PARAMS$experiment, PARAMS$scope, PARAMS$job_id,
                PARAMS$project_id, PARAMS$cdan_cv))


# ─── Threading ───────────────────────────────────────────────────────────────
# Machine-agnostic. Matches HPC pattern from i-spi/research/hpc_bayes_batch.R.
#
# Parallelism strategy:
#   - Each antigen fit uses 4 chains × 1 core/chain = 4 cores
#   - Multiple antigens run in parallel via mclapply
#   - N_PARALLEL_ANTIGENS = floor(available_cores / 4)
#   - No within-chain threading (STAN_NUM_THREADS=1) — rstan's reduce_sum
#     via TBB doesn't reliably help with typical immunoassay dataset sizes
#
# Examples:
#    4 cores →  1 antigen  × 4 chains =  4 cores (sequential)
#    8 cores →  2 antigens × 4 chains =  8 cores
#   16 cores →  4 antigens × 4 chains = 16 cores
#   32 cores →  8 antigens × 4 chains = 32 cores
#   80 cores → 20 antigens × 4 chains = 80 cores

# Detect available cores (nproc is reliable in containers)
AVAILABLE_CORES <- tryCatch(
  as.integer(system("nproc", intern = TRUE)),
  error = function(e) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) cores <- parallel::detectCores()
    if (is.na(cores)) cores <- 4L
    cores
  })

N_CHAINS          <- 4L
THREADS_PER_CHAIN <- 1L
N_ITER            <- 1000L
FAMILIES          <- c("4pl", "5pl", "gompertz")

# How many antigens to run in parallel via mclapply
# Each antigen uses N_CHAINS cores, so max parallel = floor(avail / N_CHAINS)
N_PARALLEL_ANTIGENS <- max(1L, floor(AVAILABLE_CORES / N_CHAINS))

# Lock threading: 1 thread per chain, no TBB, no OpenMP
Sys.setenv(STAN_NUM_THREADS = "1")
Sys.setenv(MC_CORES = as.character(N_CHAINS))
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(BLIS_NUM_THREADS = "1")
options(mc.cores = N_CHAINS, warn = 1)

message(sprintf("Threading: %d cores, %d parallel antigens x %d chains = %d cores used",
                AVAILABLE_CORES, N_PARALLEL_ANTIGENS, N_CHAINS,
                N_PARALLEL_ANTIGENS * N_CHAINS))


# ─── Progress Reporting ──────────────────────────────────────────────────────

PROGRESS_FILE <- file.path(PARAMS$progress_dir, sprintf("progress_%s.json", PARAMS$job_id))

write_progress <- function(total_combos, completed_combos,
                           current_experiment = "",
                           current_antigens = "",
                           experiments_done = 0L,
                           experiments_total = 0L) {
  progress <- list(
    total_combos       = total_combos,
    completed_combos   = completed_combos,
    current_experiment = current_experiment,
    current_antigens   = current_antigens,
    experiments_done   = experiments_done,
    experiments_total  = experiments_total
  )
  tryCatch(
    write(toJSON(progress, auto_unbox = TRUE), file = PROGRESS_FILE),
    error = function(e) message("Warning: could not write progress file: ", e$message)
  )
}


# ─── Inlined correct_prozone() ──────────────────────────────────────────────

correct_prozone <- function(stdframe = NULL, prop_diff = NULL, dil_scale = NULL,
                            response_variable = "mfi",
                            independent_variable = "concentration",
                            verbose = TRUE) {
  response_variable <- unique(response_variable)
  stdframe <- stdframe[!is.na(stdframe[[response_variable]]) &
                       !is.na(stdframe[[independent_variable]]), ]

  max_response <- max(stdframe[[response_variable]], na.rm = TRUE)
  logc_at_max_response <- max(
    stdframe[stdframe[[response_variable]] == max_response, ][[independent_variable]]
  )
  if (verbose) cat("Peak MFI =", max_response, "at concentration =", logc_at_max_response, "\n")

  post_peak <- stdframe[[independent_variable]] > logc_at_max_response
  if (verbose) cat("Number of points beyond the peak:", sum(post_peak), "\n")

  stdframe[post_peak, ][[response_variable]] <- max_response +
    (
      (max_response - stdframe[post_peak, ][[response_variable]]) * prop_diff /
        ((stdframe[post_peak, ][[independent_variable]] - logc_at_max_response) * dil_scale)
    )
  return(stdframe)
}


`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[[1]])) a else b


# ─── DB Connection ───────────────────────────────────────────────────────────

open_conn <- function() {
  DBI::dbConnect(
    RPostgres::Postgres(),
    dbname   = Sys.getenv("DB_NAME",    "local_madi_ispi"),
    host     = Sys.getenv("DB_HOST",    "localhost"),
    port     = as.integer(Sys.getenv("DB_PORT", "5432")),
    user     = Sys.getenv("DB_USER",    "hardik"),
    password = Sys.getenv("DB_PASSWORD", ""),
    sslmode  = Sys.getenv("DB_SSLMODE", "disable"),
    options  = "-c search_path=madi_results",
    connect_timeout = 10L
  )
}


# ─── Combo Discovery ────────────────────────────────────────────────────────

discover_combos <- function(conn, study, project_id, experiment = NULL, antigen = NULL, source_val = NULL) {
  base_sql <- "
    SELECT DISTINCT
      h.project_id,
      s.study_accession,
      s.experiment_accession,
      s.antigen,
      s.feature,
      s.source,
      COALESCE(s.wavelength, '__none__') AS wavelength
    FROM madi_results.xmap_standard s
    INNER JOIN madi_results.xmap_header h
      ON  h.study_accession      = s.study_accession
      AND h.experiment_accession = s.experiment_accession
      AND TRIM(h.plate_id)       = TRIM(s.plate_id)
    WHERE s.study_accession = $1
      AND h.project_id = $2"

  params <- list(study, project_id)
  np <- 3L

  if (!is.null(experiment) && nzchar(experiment)) {
    base_sql <- paste0(base_sql, sprintf(" AND s.experiment_accession = $%d", np))
    params <- c(params, list(experiment))
    np <- np + 1L
  }

  if (!is.null(antigen) && nzchar(antigen)) {
    base_sql <- paste0(base_sql, sprintf(" AND s.antigen = $%d", np))
    params <- c(params, list(antigen))
    np <- np + 1L
  }

  if (!is.null(source_val) && nzchar(source_val)) {
    base_sql <- paste0(base_sql, sprintf(" AND s.source = $%d", np))
    params <- c(params, list(source_val))
  }

  base_sql <- paste0(base_sql, " ORDER BY s.experiment_accession, s.antigen")
  DBI::dbGetQuery(conn, base_sql, params = params)
}


# ─── Data Fetchers ──────────────────────────────────────────────────────────

fetch_standards <- function(conn, study, experiment, antigen, wavelength, source_val = NULL) {
  base_sql <- "
    SELECT s.antigen, h.plateid, h.plate, h.nominal_sample_dilution,
           h.sample_dilution_factor, s.antibody_mfi AS mfi,
           s.dilution AS dilution_factor, s.feature, s.source,
           COALESCE(s.wavelength, '__none__') AS wavelength, h.project_id
    FROM madi_results.xmap_standard s
    INNER JOIN madi_results.xmap_header h
      ON h.study_accession = s.study_accession
      AND h.experiment_accession = s.experiment_accession
      AND TRIM(h.plate_id) = TRIM(s.plate_id)
    WHERE s.study_accession = $1 AND s.experiment_accession = $2 AND s.antigen = $3"
  params <- list(study, experiment, antigen)
  np <- 4L
  if (wavelength != "__none__") {
    base_sql <- paste0(base_sql, sprintf(" AND s.wavelength = $%d", np))
    params <- c(params, list(wavelength)); np <- np + 1L
  }
  if (!is.null(source_val) && nzchar(source_val)) {
    base_sql <- paste0(base_sql, sprintf(" AND s.source = $%d", np))
    params <- c(params, list(source_val))
  }
  DBI::dbGetQuery(conn, base_sql, params = params)
}


fetch_samples <- function(conn, study, experiment, antigen) {
  DBI::dbGetQuery(conn, "
    SELECT s.antigen, h.plateid, h.plate, h.nominal_sample_dilution,
           s.sampleid, s.antibody_mfi AS mfi, s.dilution, s.timeperiod,
           s.patientid, s.well, s.feature, s.source,
           COALESCE(s.wavelength, '__none__') AS wavelength, h.project_id
    FROM madi_results.xmap_sample s
    INNER JOIN madi_results.xmap_header h
      ON h.study_accession = s.study_accession
      AND h.experiment_accession = s.experiment_accession
      AND TRIM(h.plate_id) = TRIM(s.plate_id)
    WHERE s.study_accession = $1 AND s.experiment_accession = $2 AND s.antigen = $3",
    params = list(study, experiment, antigen))
}


fetch_blanks <- function(conn, study, experiment, antigen, wavelength, source_val = NULL) {
  base_sql <- "
    SELECT b.antigen, h.plateid, b.antibody_mfi AS mfi
    FROM madi_results.xmap_buffer b
    INNER JOIN madi_results.xmap_header h
      ON h.study_accession = b.study_accession
      AND h.experiment_accession = b.experiment_accession
      AND TRIM(h.plate_id) = TRIM(b.plate_id)
    WHERE b.study_accession = $1 AND b.experiment_accession = $2 AND b.antigen = $3
      AND UPPER(b.stype) = 'B' AND b.antibody_mfi > 0"
  params <- list(study, experiment, antigen)
  np <- 4L
  if (wavelength != "__none__") {
    base_sql <- paste0(base_sql, sprintf(" AND b.wavelength = $%d", np))
    params <- c(params, list(wavelength)); np <- np + 1L
  }
  if (!is.null(source_val) && nzchar(source_val)) {
    base_sql <- paste0(base_sql, sprintf(" AND b.source = $%d", np))
    params <- c(params, list(source_val))
  }
  DBI::dbGetQuery(conn, base_sql, params = params)
}


# ─── Helpers ─────────────────────────────────────────────────────────────────

safe_get <- function(lst, key) {
  if (is.null(lst) || is.null(lst[[key]]) || length(lst[[key]]) == 0) return(NA_real_)
  val <- lst[[key]]
  if (is.na(val)) NA_real_ else val
}


eval_forward <- function(x, a, b, c, d, g = NA, family = "5pl") {
  if (is.null(x) || is.na(x) || x <= 0) return(NA_real_)
  lx <- log(x); lc <- log(c)
  if (family %in% c("4pl", "5pl")) {
    z <- exp(b * (lx - lc)); d + (a - d) / (1 + z)^g
  } else if (family == "gompertz") {
    u <- b * (lx - lc); d + (a - d) * exp(-exp(u))
  } else NA_real_
}


compute_dydx_inflect <- function(a, b, c, d, g, x_infl, family) {
  if (is.na(x_infl) || x_infl <= 0) return(NA_real_)
  if (family %in% c("4pl", "5pl")) {
    z <- if (family == "5pl" && !is.na(g) && g > 0) 1 / g else 1
    dydlogx <- -(a - d) * (if (is.na(g)) 1 else g) * b * z / (1 + z)^((if (is.na(g)) 1 else g) + 1)
    dydlogx / x_infl
  } else if (family == "gompertz") {
    -(a - d) * b * exp(-1) / x_infl
  } else NA_real_
}


extract_plate_params_v2 <- function(assay) {
  ens <- assay$ensemble; pm <- assay$fit_obj$plate_map
  rows <- lapply(seq_len(nrow(pm)), function(p) {
    pid <- as.character(pm$plateid[p])
    fam <- ens$plate_best_family[pid]
    fit_f <- ens$fits[[fam]]; samp <- rstan::extract(fit_f$fit)
    pi <- fit_f$plate_map$plate_idx[as.character(fit_f$plate_map$plateid) == pid]
    data.frame(plateid = pid, curve_family = fam,
               a = median(samp$a[, pi]), b = median(samp$b[, pi]),
               c = exp(median(samp$log_c[, pi])), d = median(samp$d[, pi]),
               g = if ("g" %in% names(samp)) median(samp$g[, pi]) else NA_real_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


# Derive LOQ from CDAN profile at a custom CV threshold
# Mirrors the internal find_loq() in stanassay::cdan.R
find_loq_custom <- function(cdan_result, threshold) {
  if (is.null(cdan_result) || is.null(cdan_result$profile)) {
    return(list(lloq = NA_real_, uloq = NA_real_))
  }
  prof <- cdan_result$profile
  valid <- which(!is.na(prof$smoothed_cv) & prof$smoothed_cv < threshold)
  if (length(valid) == 0) return(list(lloq = NA_real_, uloq = NA_real_))
  list(
    lloq = prof$concentration[min(valid)],
    uloq = prof$concentration[max(valid)]
  )
}


compute_se_pcov_v2 <- function(plate_fit, plate_id, mfi_val) {
  NA_r <- list(raw_predicted_concentration = NA_real_, se_concentration = NA_real_,
               pcov = NA_real_, conc_lower = NA_real_, conc_upper = NA_real_)
  pred <- tryCatch(stanassay:::predict_concentration_bayesian(plate_fit, plate_id, mfi_val),
                   error = function(e) list(est = NA_real_, lower = NA_real_, upper = NA_real_))
  if (is.na(pred$est)) return(NA_r)
  pm <- plate_fit$plate_map
  pi <- pm$plate_idx[as.character(pm$plateid) == as.character(plate_id)]
  if (length(pi) == 0) return(NA_r)
  params <- tryCatch(stanassay:::.extract_curve_params(plate_fit, pi), error = function(e) NULL)
  if (is.null(params)) return(list(raw_predicted_concentration = pred$est,
                                    se_concentration = NA_real_, pcov = NA_real_,
                                    conc_lower = pred$lower, conc_upper = pred$upper))
  inv <- tryCatch(stanassay:::.curve_inverse(params, mfi_val), error = function(e) NULL)
  if (is.null(inv) || sum(inv$valid) < length(params$a) * 0.1)
    return(list(raw_predicted_concentration = pred$est, se_concentration = NA_real_,
                pcov = NA_real_, conc_lower = pred$lower, conc_upper = pred$upper))
  x_pred <- exp(inv$log_x[inv$valid])
  med <- median(x_pred); sdev <- sd(x_pred)
  list(raw_predicted_concentration = med, se_concentration = sdev,
       pcov = if (!is.na(med) && med > 0) sdev / med else NA_real_,
       conc_lower = pred$lower, conc_upper = pred$upper)
}


# ─── Fit One Antigen Combo ──────────────────────────────────────────────────

fit_one <- function(conn, study, experiment, antigen,
                    feature, source_val, wavelength, project_id,
                    cdan_cv = 20) {

  message(sprintf("  Fitting %s | %s | %s", study, experiment, antigen))

  stds_raw <- fetch_standards(conn, study, experiment, antigen, wavelength, source_val)
  if (nrow(stds_raw) == 0L) { message("  SKIP: no standards"); return(NULL) }
  samps_raw  <- fetch_samples(conn, study, experiment, antigen)
  blanks_raw <- fetch_blanks(conn, study, experiment, antigen, wavelength, source_val)

  # Look up standard_curve_concentration from xmap_antigen_family
  # This is the undiluted concentration (e.g. 10000 for tt), used as the numerator
  # in concentration = base_num / dilution_factor (same as i-spi std_curver_ui.R line 2620)
  sc_conc <- tryCatch({
    res <- DBI::dbGetQuery(conn, sprintf(
      "SELECT DISTINCT standard_curve_concentration
       FROM madi_results.xmap_antigen_family
       WHERE study_accession = '%s' AND experiment_accession = '%s' AND antigen = '%s'
         AND standard_curve_concentration IS NOT NULL
       LIMIT 1", study, experiment, antigen))
    if (nrow(res) > 0) as.numeric(res$standard_curve_concentration[1]) else NULL
  }, error = function(e) NULL)

  stds <- stds_raw |>
    mutate(plateid = as.character(plateid), mfi = as.numeric(mfi),
           dilution_factor = as.numeric(dilution_factor),
           nom = suppressWarnings(as.numeric(nominal_sample_dilution))) |>
    filter(!is.na(plateid), nchar(plateid) > 0,
           mfi > 0, !is.na(dilution_factor), dilution_factor > 0) |>
    mutate(base_num = if (!is.null(sc_conc)) sc_conc else coalesce(nom, 100000),
           concentration = base_num / dilution_factor) |>
    group_by(plateid, concentration) |>
    summarise(mfi = median(mfi, na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(mfi))

  if (!is.null(sc_conc)) {
    message(sprintf("  Using standard_curve_concentration = %g from xmap_antigen_family", sc_conc))
  } else {
    message("  No standard_curve_concentration found, falling back to nominal_sample_dilution")
  }

  if (nrow(stds) < 4L) { message("  SKIP: <4 std points"); return(NULL) }

  stds <- stds |>
    group_by(plateid) |>
    group_modify(~ correct_prozone(.x, prop_diff = 0.1, dil_scale = 2,
                                    response_variable = "mfi",
                                    independent_variable = "concentration",
                                    verbose = FALSE)) |>
    ungroup()

  samps <- samps_raw |>
    mutate(plateid = as.character(plateid), mfi = as.numeric(mfi)) |>
    filter(!is.na(plateid), mfi > 0, !is.na(mfi))
  blanks_df <- blanks_raw |>
    mutate(plateid = as.character(plateid), mfi = as.numeric(mfi)) |>
    filter(!is.na(plateid), nchar(plateid) > 0, !is.na(mfi), mfi > 0)

  assay_type <- if (wavelength != "__none__") "elisa" else "xmap"
  assay <- tryCatch({
    a <- StanAssay$new(
      std_data = stds, concentration_col = "concentration",
      response_col = "mfi", plate_col = "plateid", assay_type = assay_type,
      blank_data = if (nrow(blanks_df) > 0L) blanks_df else NULL,
      blank_response_col = "mfi", blank_plate_col = "plateid",
      apply_prozone = FALSE
    )
    suppressWarnings(a$fit_ensemble(
      families = FAMILIES, error_model = "exact", n_iter = N_ITER,
      n_chains = N_CHAINS, cores = N_CHAINS,
      threads_per_chain = THREADS_PER_CHAIN, grainsize = 1L, refresh = 0L
    ))
    a
  }, error = function(e) { message("  ERROR: ", e$message); NULL })
  if (is.null(assay)) return(NULL)

  ens <- assay$ensemble
  plate_best <- ens$plate_best_family; plate_elpd <- ens$plate_elpd
  stacking_wts <- ens$stacking_weights; global_best <- ens$best_family
  loo_comp <- ens$loo_comparison
  target_plates <- sort(unique(stds$plateid))
  params <- extract_plate_params_v2(assay)

  plate_meta <- stds_raw |>
    mutate(plateid = as.character(plateid)) |>
    select(plateid, plate, nominal_sample_dilution, sample_dilution_factor) |>
    distinct() |> group_by(plateid) |> slice(1) |> ungroup()

  # ── df_curves ──────────────────────────────────────────────────────────
  curve_rows <- lapply(target_plates, function(plt) {
    p <- params[params$plateid == plt, ]
    if (nrow(p) == 0) return(NULL)
    fam <- p$curve_family; meta <- plate_meta[plate_meta$plateid == plt, ]

    cdan     <- tryCatch(assay$compute_cdan(plt), error = function(e) NULL)
    d2       <- tryCatch(assay$compute_second_derivative(plt, n_grid = 200L), error = function(e) NULL)
    lod_res  <- tryCatch(assay$compute_lod(plt), error = function(e) NULL)
    lrdl_res <- tryCatch(assay$compute_lrdl(plt), error = function(e) NULL)
    uod_res  <- tryCatch(assay$compute_uod(plt), error = function(e) NULL)
    urdl_res <- tryCatch(assay$compute_urdl(plt), error = function(e) NULL)
    infl     <- tryCatch(assay$compute_inflection_point(plt), error = function(e) NULL)

    inflect_x <- safe_get(infl, "x_median"); inflect_y <- safe_get(infl, "y_median")
    dydx <- compute_dydx_inflect(p$a, p$b, p$c, p$d, p$g, inflect_x, fam)
    fwd <- function(x) eval_forward(x, p$a, p$b, p$c, p$d, p$g, fam)

    lloq_20 <- safe_get(cdan, "lloq_20"); uloq_20 <- safe_get(cdan, "uloq_20")
    lloq_15 <- safe_get(cdan, "lloq_15"); uloq_15 <- safe_get(cdan, "uloq_15")
    custom_loq <- find_loq_custom(cdan, cdan_cv)
    lloq_custom <- custom_loq$lloq; uloq_custom <- custom_loq$uloq
    lo2d <- safe_get(d2, "lo2d"); uo2d <- safe_get(d2, "uo2d")
    lod_conc <- safe_get(lod_res, "lod"); lod_y <- safe_get(lod_res, "threshold_mfi_median")
    lrdl_conc <- safe_get(lrdl_res, "lrdl")
    uod_conc <- safe_get(uod_res, "uod"); urdl_conc <- safe_get(urdl_res, "urdl")
    elpd_row <- if (plt %in% rownames(plate_elpd)) plate_elpd[plt, ] else rep(NA_real_, length(FAMILIES))

    data.frame(
      project_id = project_id, study_accession = study, experiment_accession = experiment,
      plateid = plt,
      plate = if (nrow(meta) > 0) as.character(meta$plate[1]) else NA_character_,
      nominal_sample_dilution = if (nrow(meta) > 0) as.character(meta$nominal_sample_dilution[1]) else NA_character_,
      sample_dilution_factor = if (nrow(meta) > 0) as.numeric(meta$sample_dilution_factor[1]) else NA_real_,
      feature = feature, antigen = antigen, source = source_val, wavelength = wavelength,
      apply_prozone = FALSE, curve_family = fam,
      a = p$a, b = p$b, c = p$c, d = p$d, g = p$g,
      lloq = lloq_20, uloq = uloq_20, lloq_y = fwd(lloq_20), uloq_y = fwd(uloq_20),
      lloq_15 = lloq_15, uloq_15 = uloq_15,
      lloq_custom = lloq_custom, uloq_custom = uloq_custom,
      cdan_cv_threshold = cdan_cv,
      lo2d = lo2d, uo2d = uo2d, lo2d_y = fwd(lo2d), uo2d_y = fwd(uo2d),
      lod = lod_conc, lod_y = lod_y, lrdl = lrdl_conc,
      uod = uod_conc, uod_y = fwd(uod_conc), urdl = urdl_conc,
      inflect_x = inflect_x, inflect_y = inflect_y,
      inflect_x_lower = safe_get(infl, "x_lower"), inflect_x_upper = safe_get(infl, "x_upper"),
      dydx_inflect = dydx, plate_best_family = fam,
      plate_elpd_4pl = elpd_row["4pl"], plate_elpd_5pl = elpd_row["5pl"],
      plate_elpd_gompertz = elpd_row["gompertz"],
      global_stacking_4pl = stacking_wts["4pl"], global_stacking_5pl = stacking_wts["5pl"],
      global_stacking_gompertz = stacking_wts["gompertz"],
      global_best_family = global_best,
      cdan_method = if (!is.null(cdan) && !is.null(cdan$method)) cdan$method else NA_character_,
      stringsAsFactors = FALSE, row.names = NULL)
  })
  df_curves_local <- do.call(rbind, Filter(Negate(is.null), curve_rows))

  # ── df_samples ─────────────────────────────────────────────────────────
  df_samples_local <- NULL
  if (nrow(samps) > 0L) {
    message(sprintf("  Back-calculating %d samples...", nrow(samps)))
    se_rows <- lapply(seq_len(nrow(samps)), function(i) {
      plt <- samps$plateid[i]
      pf <- tryCatch(assay$.__enclos_env__$private$.get_plate_fit_obj(plt),
                     error = function(e) assay$fit_obj)
      compute_se_pcov_v2(pf, plt, samps$mfi[i])
    })
    df_samples_local <- samps |>
      mutate(project_id = project_id, study_accession = study, experiment_accession = experiment,
             feature = feature, source = source_val, wavelength = wavelength,
             raw_predicted_concentration = vapply(se_rows, `[[`, numeric(1), "raw_predicted_concentration"),
             se_concentration = vapply(se_rows, `[[`, numeric(1), "se_concentration"),
             pcov = vapply(se_rows, `[[`, numeric(1), "pcov"),
             conc_lower = vapply(se_rows, `[[`, numeric(1), "conc_lower"),
             conc_upper = vapply(se_rows, `[[`, numeric(1), "conc_upper"))
    df_samples_local$gate_class <- NA_character_
    for (plt in target_plates) {
      cn <- tryCatch(assay$compute_cdan(plt), error = function(e) NULL)
      if (is.null(cn)) next
      ll <- safe_get(cn, "lloq_20"); ul <- safe_get(cn, "uloq_20")
      if (is.na(ll) || is.na(ul)) next
      mask <- df_samples_local$plateid == plt & !is.na(df_samples_local$raw_predicted_concentration)
      cc <- df_samples_local$raw_predicted_concentration[mask]
      df_samples_local$gate_class[mask] <- ifelse(cc < ll, "Below LLOQ",
                                                   ifelse(cc > ul, "Above ULOQ", "Within Range"))
    }
    df_samples_local <- df_samples_local |>
      select(any_of(c("project_id", "study_accession", "experiment_accession",
                       "plateid", "plate", "antigen", "feature",
                       "nominal_sample_dilution", "source", "wavelength",
                       "timeperiod", "patientid", "well", "sampleid", "dilution",
                       "mfi", "raw_predicted_concentration", "se_concentration", "pcov",
                       "conc_lower", "conc_upper", "gate_class")))
  }

  # ── df_ensemble (with 95% CIs for each parameter) ───────────────────
  fams <- names(ens$stacking_weights); erows <- list()
  for (fn in fams) {
    fit_f <- ens$fits[[fn]]; samp <- rstan::extract(fit_f$fit)
    for (pid in target_plates) {
      pi <- fit_f$plate_map$plate_idx[as.character(fit_f$plate_map$plateid) == pid]
      if (length(pi) == 0) next
      hg <- "g" %in% names(samp)

      a_draws <- samp$a[, pi]
      b_draws <- samp$b[, pi]
      c_draws <- exp(samp$log_c[, pi])
      d_draws <- samp$d[, pi]
      g_draws <- if (hg) samp$g[, pi] else rep(NA_real_, length(a_draws))

      erows[[length(erows) + 1]] <- data.frame(
        project_id = project_id, study_accession = study, experiment_accession = experiment,
        plateid = pid, antigen = antigen, family = fn,
        plate_elpd = if (pid %in% rownames(plate_elpd)) plate_elpd[pid, fn] else NA_real_,
        is_plate_best = identical(plate_best[pid], fn),
        global_stacking_weight = stacking_wts[fn], is_global_best = identical(global_best, fn),
        a       = median(a_draws),
        a_lower = quantile(a_draws, 0.025),
        a_upper = quantile(a_draws, 0.975),
        b       = median(b_draws),
        b_lower = quantile(b_draws, 0.025),
        b_upper = quantile(b_draws, 0.975),
        c       = median(c_draws),
        c_lower = quantile(c_draws, 0.025),
        c_upper = quantile(c_draws, 0.975),
        d       = median(d_draws),
        d_lower = quantile(d_draws, 0.025),
        d_upper = quantile(d_draws, 0.975),
        g       = if (hg) median(g_draws) else NA_real_,
        g_lower = if (hg) quantile(g_draws, 0.025) else NA_real_,
        g_upper = if (hg) quantile(g_draws, 0.975) else NA_real_,
        stringsAsFactors = FALSE, row.names = NULL)
    }
  }
  df_ensemble_local <- do.call(rbind, erows)

  # ── plots: generated interactive plotly objects for each plate ──────
  plots_local <- list()
  for (plt in target_plates) {
    p_obj <- tryCatch(
      assay$plot(plt, sample_data = df_samples_local),
      error = function(e) {
        warning(sprintf("Plot failed for %s: %s", plt, e$message))
        NULL
      }
    )
    if (!is.null(p_obj)) plots_local[[plt]] <- p_obj
  }

  # ── curve_grid + cdan_grid: precomputed plot data per plate ────────
  curve_grid_rows <- list()
  cdan_grid_rows  <- list()
  for (plt in target_plates) {
    p <- params[params$plateid == plt, ]
    if (nrow(p) == 0) next
    meta <- plate_meta[plate_meta$plateid == plt, ]

    # Curve grid: posterior median + 95% CI ribbon from build_plot_data()
    pd <- tryCatch(
      stanassay::build_plot_data(assay$fit_obj, plt, stds, n_points = 100L, n_draws = 200L),
      error = function(e) { warning(sprintf("build_plot_data failed for %s: %s", plt, e$message)); NULL }
    )
    if (!is.null(pd)) {
      curve_grid_rows[[length(curve_grid_rows) + 1]] <- data.frame(
        project_id = project_id, study_accession = study, experiment_accession = experiment,
        plateid = plt, antigen = antigen, source = source_val,
        log10_conc    = log10(pd$curve_df$concentration),
        concentration = pd$curve_df$concentration,
        mfi_median    = pd$curve_df$mfi_median,
        mfi_lower_95  = pd$curve_df$mfi_lower_95,
        mfi_upper_95  = pd$curve_df$mfi_upper_95,
        stringsAsFactors = FALSE
      )
    }

    # CDAN grid: precision profile from compute_cdan
    cdan_res <- tryCatch(
      stanassay::compute_cdan_precision_profile(assay$fit_obj, plt, n_grid = 200L),
      error = function(e) NULL
    )
    if (!is.null(cdan_res) && !is.null(cdan_res$profile)) {
      prof <- cdan_res$profile
      cdan_grid_rows[[length(cdan_grid_rows) + 1]] <- data.frame(
        project_id = project_id, study_accession = study, experiment_accession = experiment,
        plateid = plt, antigen = antigen, source = source_val,
        log10_conc    = prof$log10_conc,
        concentration = prof$concentration,
        cv_percent    = prof$cv_percent,
        smoothed_cv   = prof$smoothed_cv,
        median_mu     = prof$median_mu,
        stringsAsFactors = FALSE
      )
    }
  }
  df_curve_grid <- if (length(curve_grid_rows) > 0) do.call(rbind, curve_grid_rows) else data.frame()
  df_cdan_grid  <- if (length(cdan_grid_rows) > 0) do.call(rbind, cdan_grid_rows) else data.frame()

  message(sprintf("  Done: %d curves, %d samples, %d ensemble, %d plots, %d curve_grid, %d cdan_grid",
                  nrow(df_curves_local), nrow(df_samples_local) %||% 0L,
                  nrow(df_ensemble_local), length(plots_local),
                  nrow(df_curve_grid), nrow(df_cdan_grid)))
  list(curves = df_curves_local, samples = df_samples_local,
       ensemble = df_ensemble_local, loo_comp = loo_comp,
       plots = plots_local,
       curve_grid = df_curve_grid, cdan_grid = df_cdan_grid)
}


# ─── Save to PostgreSQL ──────────────────────────────────────────────────────

upsert_on_conflict <- function(conn, schema, table, df, nk_cols) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(0L))

  # Get DB column names and filter df to only those that exist
  db_cols <- DBI::dbGetQuery(conn, sprintf(
    "SELECT column_name FROM information_schema.columns
     WHERE table_schema = '%s' AND table_name = '%s'", schema, table
  ))$column_name
  df <- df[, names(df) %in% db_cols, drop = FALSE]

  if (ncol(df) == 0L) { message("  No matching columns for ", table); return(invisible(0L)) }

  # Write to temp table, then INSERT ... ON CONFLICT DO UPDATE
  tmp <- paste0("tmp_", table)
  DBI::dbWriteTable(conn, tmp, df, overwrite = TRUE, row.names = FALSE)

  cols <- names(df)
  col_list <- paste(sprintf('"%s"', cols), collapse = ", ")
  nk_list <- paste(sprintf('"%s"', nk_cols), collapse = ", ")
  update_cols <- setdiff(cols, nk_cols)
  update_set <- paste(
    sprintf('"%s" = EXCLUDED."%s"', update_cols, update_cols),
    collapse = ", "
  )

  sql <- sprintf(
    'INSERT INTO %s.%s (%s) SELECT %s FROM %s ON CONFLICT (%s) DO UPDATE SET %s',
    schema, table, col_list, col_list, tmp, nk_list, update_set
  )
  n <- DBI::dbExecute(conn, sql)
  DBI::dbRemoveTable(conn, tmp)
  message(sprintf("  Upserted %d rows into %s.%s", n, schema, table))
  invisible(n)
}


save_to_db <- function(df_curves, df_samples, df_ensemble, job_id,
                       df_curve_grid = NULL, df_cdan_grid = NULL) {
  conn <- open_conn()
  on.exit(DBI::dbDisconnect(conn), add = TRUE)

  # Add job_id (created_at is handled by DB DEFAULT NOW())
  if (!is.null(df_curves) && nrow(df_curves) > 0) {
    df_curves$job_id <- job_id
  }
  if (!is.null(df_samples) && nrow(df_samples) > 0) {
    df_samples$job_id <- job_id
  }
  if (!is.null(df_ensemble) && nrow(df_ensemble) > 0) {
    df_ensemble$job_id <- job_id
  }

  curves_nk <- c("project_id", "study_accession", "experiment_accession",
                  "plateid", "plate", "nominal_sample_dilution",
                  "source", "wavelength", "antigen", "feature")
  samples_nk <- c(curves_nk, "patientid", "timeperiod", "sampleid", "dilution")
  ensemble_nk <- c("project_id", "study_accession", "experiment_accession",
                    "plateid", "antigen", "family")

  # 1. Upsert bayes_curves (parent)
  message("  Saving bayes_curves...")
  upsert_on_conflict(conn, "madi_results", "bayes_curves", df_curves, curves_nk)

  # 2. Fetch back bayes_curves_id for FK joins
  if (nrow(df_curves) > 0) {
    study <- unique(df_curves$study_accession)[1]
    proj <- unique(df_curves$project_id)[1]
    exp_list <- unique(df_curves$experiment_accession)
    exp_sql <- paste0("'", paste(exp_list, collapse = "','"), "'")

    lookup <- DBI::dbGetQuery(conn, sprintf(
      "SELECT bayes_curves_id, project_id, study_accession, experiment_accession,
              plateid, plate, nominal_sample_dilution, source, wavelength, antigen, feature
       FROM madi_results.bayes_curves
       WHERE project_id = %d AND study_accession = '%s'
         AND experiment_accession IN (%s)",
      proj, study, exp_sql
    ))

    # FK join for samples
    if (!is.null(df_samples) && nrow(df_samples) > 0) {
      df_samples <- dplyr::inner_join(df_samples, lookup, by = curves_nk)
      message("  Saving bayes_samples...")
      upsert_on_conflict(conn, "madi_results", "bayes_samples", df_samples, samples_nk)
    }

    # FK join for ensemble — ensemble NK doesn't include source/wavelength/etc,
    # so we need a deduplicated lookup (one bayes_curves_id per plateid+antigen)
    if (!is.null(df_ensemble) && nrow(df_ensemble) > 0) {
      ensemble_lookup <- lookup %>%
        dplyr::distinct(project_id, study_accession, experiment_accession,
                        plateid, antigen, .keep_all = TRUE) %>%
        dplyr::select(project_id, study_accession, experiment_accession,
                      plateid, antigen, bayes_curves_id)
      df_ensemble <- dplyr::inner_join(
        df_ensemble, ensemble_lookup,
        by = c("project_id", "study_accession", "experiment_accession", "plateid", "antigen")
      )
      message("  Saving bayes_ensemble...")
      upsert_on_conflict(conn, "madi_results", "bayes_ensemble", df_ensemble, ensemble_nk)
    }

    # Save curve_grid (delete old + insert — no NK, just FK)
    if (!is.null(df_curve_grid) && nrow(df_curve_grid) > 0) {
      df_curve_grid$job_id <- job_id
      # FK join: need bayes_curves_id per (plateid, antigen, source)
      grid_lookup <- lookup %>%
        dplyr::select(project_id, study_accession, experiment_accession,
                      plateid, antigen, source, bayes_curves_id) %>%
        dplyr::distinct()
      df_curve_grid <- dplyr::inner_join(
        df_curve_grid, grid_lookup,
        by = c("project_id", "study_accession", "experiment_accession",
               "plateid", "antigen", "source")
      )
      # Delete old grid rows for these plates, then insert fresh
      plateids_sql <- paste0("'", unique(df_curve_grid$plateid), "'", collapse = ",")
      sources_sql <- paste0("'", unique(df_curve_grid$source), "'", collapse = ",")
      DBI::dbExecute(conn, sprintf(
        "DELETE FROM madi_results.bayes_curve_grid
         WHERE project_id = %d AND study_accession = '%s'
           AND experiment_accession IN (%s)
           AND plateid IN (%s) AND antigen = '%s' AND source IN (%s)",
        proj, study, exp_sql, plateids_sql, unique(df_curve_grid$antigen)[1], sources_sql
      ))
      # Filter to DB columns
      db_cols <- DBI::dbGetQuery(conn,
        "SELECT column_name FROM information_schema.columns
         WHERE table_schema = 'madi_results' AND table_name = 'bayes_curve_grid'"
      )$column_name
      df_cg <- df_curve_grid[, names(df_curve_grid) %in% db_cols, drop = FALSE]
      DBI::dbWriteTable(conn, "tmp_bayes_curve_grid", df_cg, overwrite = TRUE, row.names = FALSE)
      col_list <- paste(sprintf('"%s"', names(df_cg)), collapse = ", ")
      DBI::dbExecute(conn, sprintf(
        "INSERT INTO madi_results.bayes_curve_grid (%s) SELECT %s FROM tmp_bayes_curve_grid",
        col_list, col_list
      ))
      DBI::dbRemoveTable(conn, "tmp_bayes_curve_grid")
      message(sprintf("  Saved %d rows into madi_results.bayes_curve_grid", nrow(df_cg)))
    }

    # Save cdan_grid (delete old + insert — same pattern)
    if (!is.null(df_cdan_grid) && nrow(df_cdan_grid) > 0) {
      df_cdan_grid$job_id <- job_id
      grid_lookup2 <- lookup %>%
        dplyr::select(project_id, study_accession, experiment_accession,
                      plateid, antigen, source, bayes_curves_id) %>%
        dplyr::distinct()
      df_cdan_grid <- dplyr::inner_join(
        df_cdan_grid, grid_lookup2,
        by = c("project_id", "study_accession", "experiment_accession",
               "plateid", "antigen", "source")
      )
      plateids_sql2 <- paste0("'", unique(df_cdan_grid$plateid), "'", collapse = ",")
      sources_sql2 <- paste0("'", unique(df_cdan_grid$source), "'", collapse = ",")
      DBI::dbExecute(conn, sprintf(
        "DELETE FROM madi_results.bayes_cdan_grid
         WHERE project_id = %d AND study_accession = '%s'
           AND experiment_accession IN (%s)
           AND plateid IN (%s) AND antigen = '%s' AND source IN (%s)",
        proj, study, exp_sql, plateids_sql2, unique(df_cdan_grid$antigen)[1], sources_sql2
      ))
      db_cols2 <- DBI::dbGetQuery(conn,
        "SELECT column_name FROM information_schema.columns
         WHERE table_schema = 'madi_results' AND table_name = 'bayes_cdan_grid'"
      )$column_name
      df_cd <- df_cdan_grid[, names(df_cdan_grid) %in% db_cols2, drop = FALSE]
      DBI::dbWriteTable(conn, "tmp_bayes_cdan_grid", df_cd, overwrite = TRUE, row.names = FALSE)
      col_list2 <- paste(sprintf('"%s"', names(df_cd)), collapse = ", ")
      DBI::dbExecute(conn, sprintf(
        "INSERT INTO madi_results.bayes_cdan_grid (%s) SELECT %s FROM tmp_bayes_cdan_grid",
        col_list2, col_list2
      ))
      DBI::dbRemoveTable(conn, "tmp_bayes_cdan_grid")
      message(sprintf("  Saved %d rows into madi_results.bayes_cdan_grid", nrow(df_cd)))
    }
  }

  message("  DB save complete.")
}


# =============================================================================
# MAIN — Parallel antigens via mclapply, save to DB after each batch
#
# Same pattern as i-spi/research/hpc_bayes_batch.R:
#   - mclapply runs N_PARALLEL_ANTIGENS at once
#   - Each antigen uses 4 cores (4 chains × 1 thread)
#   - After each parallel batch completes, save all results to DB
#   - Progress updated per batch
# =============================================================================

conn <- open_conn()
message("Connected to database")

combos <- discover_combos(
  conn,
  study      = PARAMS$study,
  project_id = PARAMS$project_id,
  experiment = if (nzchar(PARAMS$experiment)) PARAMS$experiment else NULL,
  antigen    = if (nzchar(PARAMS$antigen)) PARAMS$antigen else NULL,
  source_val = if (nzchar(PARAMS$source)) PARAMS$source else NULL
)
DBI::dbDisconnect(conn)

if (nrow(combos) == 0L) {
  message("No combinations found. Exiting.")
  quit(status = 0)
}

experiments <- unique(combos$experiment_accession)
total_combos <- nrow(combos)
completed_combos <- 0L
experiments_done <- 0L
experiments_total <- length(experiments)

message(sprintf("\n==== Batch Worker ===="))
message(sprintf("  Study: %s  |  Project: %d", PARAMS$study, PARAMS$project_id))
message(sprintf("  Scope: %s", PARAMS$scope))
message(sprintf("  %d combos across %d experiments", total_combos, experiments_total))
message(sprintf("  Chains: %d  |  Parallel antigens: %d  |  CDAN CV: %.1f%%\n",
                N_CHAINS, N_PARALLEL_ANTIGENS, PARAMS$cdan_cv))

write_progress(total_combos, 0L, "", "", experiments_done, experiments_total)

batch_start <- Sys.time()

for (exp_name in experiments) {
  exp_combos <- combos[combos$experiment_accession == exp_name, ]
  n_ag <- nrow(exp_combos)
  message(sprintf("\n======== %s (%d antigens, %d parallel) ========",
                  exp_name, n_ag, min(n_ag, N_PARALLEL_ANTIGENS)))

  antigen_labels <- paste(exp_combos$antigen, collapse = ", ")
  write_progress(total_combos, completed_combos, exp_name, antigen_labels,
                 experiments_done, experiments_total)

  # Run antigens in parallel via mclapply (same as hpc_bayes_batch.R)
  t0 <- Sys.time()
  results_list <- mclapply(seq_len(n_ag), function(i) {
    r <- exp_combos[i, ]
    conn_local <- open_conn()
    on.exit(DBI::dbDisconnect(conn_local))
    t_start <- Sys.time()
    result <- tryCatch(
      fit_one(conn_local, r$study_accession, r$experiment_accession, r$antigen,
              r$feature, r$source, r$wavelength, r$project_id,
              cdan_cv = PARAMS$cdan_cv),
      error = function(e) { message("  FATAL: ", e$message); NULL }
    )
    el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    if (!is.null(result)) {
      message(sprintf("  OK %s | %s %.1fs (%d curves, %d samples)",
                      exp_name, r$antigen, el,
                      nrow(result$curves), nrow(result$samples) %||% 0L))
    } else {
      message(sprintf("  FAIL %s | %s %.1fs", exp_name, r$antigen, el))
    }
    result
  }, mc.cores = N_PARALLEL_ANTIGENS)
  batch_el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Save all results from this batch to DB
  n_ok <- 0L
  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    if (is.null(result)) next
    r <- exp_combos[i, ]
    combo_label <- sprintf("%s | %s", exp_name, r$antigen)
    tryCatch({
      message(sprintf("  Saving %s to DB...", combo_label))
      save_to_db(result$curves, result$samples, result$ensemble, PARAMS$job_id,
                 df_curve_grid = result$curve_grid, df_cdan_grid = result$cdan_grid)
      n_ok <- n_ok + 1L
    }, error = function(e) {
      message(sprintf("  DB SAVE ERROR for %s: %s", combo_label, e$message))
    })
  }

  completed_combos <- completed_combos + n_ag
  message(sprintf("  Batch done: %d/%d saved in %.1fs", n_ok, n_ag, batch_el))

  experiments_done <- experiments_done + 1L
  write_progress(total_combos, completed_combos, exp_name, "done",
                 experiments_done, experiments_total)
}

batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))

# Final progress
write_progress(total_combos, total_combos, "", "",
               experiments_total, experiments_total)

message(sprintf("\n==== BATCH COMPLETE — %d/%d combos in %.1f minutes ====",
                completed_combos, total_combos, batch_elapsed))
