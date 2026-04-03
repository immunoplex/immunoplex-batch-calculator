-- =============================================================================
-- 002_create_bayes_plot_data.sql
--
-- Stores pre-computed plot grids so the Bayesian curve plot can be fully
-- reconstructed from DB without needing posterior draws.
--
-- Two tables:
--   bayes_curve_grid  — fitted curve + 95% CI ribbon (100 points per plate)
--   bayes_cdan_grid   — CDAN precision profile (200 points per plate)
--
-- Run: psql -d local_madi_ispi -f migrations/002_create_bayes_plot_data.sql
-- =============================================================================

SET search_path TO madi_results, public;

-- ─── bayes_curve_grid ─────────────────────────────────────────────────────
-- Pre-computed curve: posterior median + 95% CI at log-spaced x-grid.
-- Produced by build_plot_data() in stanassay — 100 points per plate.

CREATE TABLE IF NOT EXISTS madi_results.bayes_curve_grid (
    bayes_curve_grid_id   BIGSERIAL PRIMARY KEY,

    project_id            INTEGER        NOT NULL DEFAULT 0,
    study_accession       VARCHAR(15)    NOT NULL,
    experiment_accession  VARCHAR(15)    NOT NULL,
    plateid               VARCHAR(100)   NOT NULL,
    antigen               VARCHAR(64)    NOT NULL DEFAULT '__none__',
    source                VARCHAR(25)    NOT NULL DEFAULT '__none__',

    -- Grid point
    log10_conc            NUMERIC        NOT NULL,
    concentration         NUMERIC        NOT NULL,
    mfi_median            NUMERIC        NOT NULL,
    mfi_lower_95          NUMERIC,
    mfi_upper_95          NUMERIC,

    bayes_curves_id       BIGINT REFERENCES madi_results.bayes_curves(bayes_curves_id),
    job_id                VARCHAR(64),
    created_at            TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_bayes_curve_grid_plate
    ON madi_results.bayes_curve_grid (plateid, antigen, source);


-- ─── bayes_cdan_grid ──────────────────────────────────────────────────────
-- CDAN precision profile: smoothed CV% vs concentration.
-- Produced by compute_cdan_precision_profile() — 200 points per plate.

CREATE TABLE IF NOT EXISTS madi_results.bayes_cdan_grid (
    bayes_cdan_grid_id    BIGSERIAL PRIMARY KEY,

    project_id            INTEGER        NOT NULL DEFAULT 0,
    study_accession       VARCHAR(15)    NOT NULL,
    experiment_accession  VARCHAR(15)    NOT NULL,
    plateid               VARCHAR(100)   NOT NULL,
    antigen               VARCHAR(64)    NOT NULL DEFAULT '__none__',
    source                VARCHAR(25)    NOT NULL DEFAULT '__none__',

    -- Grid point
    log10_conc            NUMERIC        NOT NULL,
    concentration         NUMERIC        NOT NULL,
    cv_percent            NUMERIC,
    smoothed_cv           NUMERIC,
    median_mu             NUMERIC,

    bayes_curves_id       BIGINT REFERENCES madi_results.bayes_curves(bayes_curves_id),
    job_id                VARCHAR(64),
    created_at            TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_bayes_cdan_grid_plate
    ON madi_results.bayes_cdan_grid (plateid, antigen, source);


-- ─── Done ─────────────────────────────────────────────────────────────────
-- Verify:
--   \dt madi_results.bayes_*_grid
