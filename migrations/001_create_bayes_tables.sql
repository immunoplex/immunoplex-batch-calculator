-- =============================================================================
-- 001_create_bayes_tables.sql
--
-- Creates the Bayesian standard curve results tables in madi_results schema.
-- Idempotent: safe to run multiple times (IF NOT EXISTS).
--
-- Tables:
--   bayes_curves   — Curve params + diagnostics (1 row per plate × antigen)
--   bayes_samples  — Sample concentrations with credible intervals
--   bayes_ensemble — Per-family model comparison + parameter posteriors
--
-- Run: psql -d local_madi_ispi -f 001_create_bayes_tables.sql
-- Prod: psql -h <host> -U <user> -d <db> -f 001_create_bayes_tables.sql
-- =============================================================================

SET search_path TO madi_results, public;

-- ─── bayes_curves ───────────────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS madi_results.bayes_curves (
    bayes_curves_id           BIGSERIAL PRIMARY KEY,

    -- Natural key (matches best_glance_all NK pattern)
    project_id                INTEGER        NOT NULL DEFAULT 0,
    study_accession           VARCHAR(15)    NOT NULL,
    experiment_accession      VARCHAR(15)    NOT NULL,
    plateid                   VARCHAR(100)   NOT NULL,
    plate                     VARCHAR(40)    NOT NULL,
    nominal_sample_dilution   VARCHAR(128)   NOT NULL,
    feature                   VARCHAR(15)    NOT NULL DEFAULT '__none__',
    antigen                   VARCHAR(64)    NOT NULL DEFAULT '__none__',
    source                    VARCHAR(25)    NOT NULL DEFAULT '__none__',
    wavelength                VARCHAR(15)    NOT NULL DEFAULT '__none__',

    -- Plate metadata
    sample_dilution_factor    NUMERIC,
    apply_prozone             BOOLEAN,

    -- Model selection (Bayesian-specific)
    curve_family              VARCHAR(20),
    plate_best_family         VARCHAR(20),
    global_best_family        VARCHAR(20),

    -- Curve parameters (posterior median)
    a                         NUMERIC,
    b                         NUMERIC,
    c                         NUMERIC,
    d                         NUMERIC,
    g                         NUMERIC,

    -- CDAN precision-based LOQ (20% CV)
    lloq                      NUMERIC,
    uloq                      NUMERIC,
    lloq_y                    NUMERIC,
    uloq_y                    NUMERIC,

    -- CDAN LOQ (15% CV)
    lloq_15                   NUMERIC,
    uloq_15                   NUMERIC,

    -- CDAN LOQ (custom CV — user-specified threshold)
    lloq_custom               NUMERIC,
    uloq_custom               NUMERIC,
    cdan_cv_threshold         NUMERIC,

    -- Second derivative limits
    lo2d                      NUMERIC,
    uo2d                      NUMERIC,
    lo2d_y                    NUMERIC,
    uo2d_y                    NUMERIC,

    -- Detection limits
    lod                       NUMERIC,
    lod_y                     NUMERIC,
    lrdl                      NUMERIC,
    uod                       NUMERIC,
    uod_y                     NUMERIC,
    urdl                      NUMERIC,

    -- Inflection point (with credible intervals)
    inflect_x                 NUMERIC,
    inflect_y                 NUMERIC,
    inflect_x_lower           NUMERIC,
    inflect_x_upper           NUMERIC,
    dydx_inflect              NUMERIC,

    -- Ensemble: per-plate ELPD
    plate_elpd_4pl            NUMERIC,
    plate_elpd_5pl            NUMERIC,
    plate_elpd_gompertz       NUMERIC,

    -- Ensemble: global stacking weights
    global_stacking_4pl       NUMERIC,
    global_stacking_5pl       NUMERIC,
    global_stacking_gompertz  NUMERIC,

    -- CDAN method used
    cdan_method               VARCHAR(50),

    -- Audit
    created_at                TIMESTAMPTZ DEFAULT NOW(),
    job_id                    VARCHAR(64),

    CONSTRAINT bayes_curves_nk UNIQUE (
        project_id, study_accession, experiment_accession,
        plateid, plate, nominal_sample_dilution,
        source, wavelength, antigen, feature
    )
);

-- Indexes for common query patterns
CREATE INDEX IF NOT EXISTS idx_bayes_curves_study
    ON madi_results.bayes_curves (study_accession, experiment_accession);
CREATE INDEX IF NOT EXISTS idx_bayes_curves_project_study
    ON madi_results.bayes_curves (project_id, study_accession);


-- ─── bayes_samples ──────────────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS madi_results.bayes_samples (
    bayes_samples_id              BIGSERIAL PRIMARY KEY,

    -- Natural key
    project_id                    INTEGER        NOT NULL DEFAULT 0,
    study_accession               VARCHAR(15)    NOT NULL,
    experiment_accession          VARCHAR(15)    NOT NULL,
    plateid                       VARCHAR(100)   NOT NULL,
    plate                         VARCHAR(40)    NOT NULL,
    nominal_sample_dilution       VARCHAR(128)   NOT NULL,
    feature                       VARCHAR(15)    NOT NULL DEFAULT '__none__',
    antigen                       VARCHAR(64)    NOT NULL DEFAULT '__none__',
    source                        VARCHAR(25)    NOT NULL DEFAULT '__none__',
    wavelength                    VARCHAR(15)    NOT NULL DEFAULT '__none__',

    -- Sample identification
    patientid                     VARCHAR(15),
    timeperiod                    VARCHAR(40),
    well                          VARCHAR(6),
    sampleid                      VARCHAR(15)    NOT NULL,
    dilution                      NUMERIC        NOT NULL,

    -- Response
    mfi                           NUMERIC,

    -- Bayesian concentration estimates
    raw_predicted_concentration   NUMERIC,
    se_concentration              NUMERIC,
    pcov                          NUMERIC,
    conc_lower                    NUMERIC,
    conc_upper                    NUMERIC,

    -- Classification
    gate_class                    VARCHAR(50),

    -- FK to parent curves table
    bayes_curves_id               BIGINT REFERENCES madi_results.bayes_curves(bayes_curves_id),

    -- Audit
    created_at                    TIMESTAMPTZ DEFAULT NOW(),
    job_id                        VARCHAR(64),

    CONSTRAINT bayes_samples_nk UNIQUE (
        project_id, study_accession, experiment_accession,
        plateid, plate, nominal_sample_dilution,
        source, wavelength, antigen, feature,
        patientid, timeperiod, sampleid, dilution
    )
);

CREATE INDEX IF NOT EXISTS idx_bayes_samples_study
    ON madi_results.bayes_samples (study_accession, experiment_accession);
CREATE INDEX IF NOT EXISTS idx_bayes_samples_curves_fk
    ON madi_results.bayes_samples (bayes_curves_id);
CREATE INDEX IF NOT EXISTS idx_bayes_samples_patient
    ON madi_results.bayes_samples (study_accession, patientid, antigen);


-- ─── bayes_ensemble ─────────────────────────────────────────────────────────

CREATE TABLE IF NOT EXISTS madi_results.bayes_ensemble (
    bayes_ensemble_id             BIGSERIAL PRIMARY KEY,

    -- Natural key
    project_id                    INTEGER        NOT NULL DEFAULT 0,
    study_accession               VARCHAR(15)    NOT NULL,
    experiment_accession          VARCHAR(15)    NOT NULL,
    plateid                       VARCHAR(100)   NOT NULL,
    antigen                       VARCHAR(64)    NOT NULL DEFAULT '__none__',
    family                        VARCHAR(20)    NOT NULL,

    -- Model comparison
    plate_elpd                    NUMERIC,
    is_plate_best                 BOOLEAN,
    global_stacking_weight        NUMERIC,
    is_global_best                BOOLEAN,

    -- Parameter posteriors (median + 95% CI)
    a                             NUMERIC,
    a_lower                       NUMERIC,
    a_upper                       NUMERIC,
    b                             NUMERIC,
    b_lower                       NUMERIC,
    b_upper                       NUMERIC,
    c                             NUMERIC,
    c_lower                       NUMERIC,
    c_upper                       NUMERIC,
    d                             NUMERIC,
    d_lower                       NUMERIC,
    d_upper                       NUMERIC,
    g                             NUMERIC,
    g_lower                       NUMERIC,
    g_upper                       NUMERIC,

    -- FK to parent curves table
    bayes_curves_id               BIGINT REFERENCES madi_results.bayes_curves(bayes_curves_id),

    -- Audit
    created_at                    TIMESTAMPTZ DEFAULT NOW(),
    job_id                        VARCHAR(64),

    CONSTRAINT bayes_ensemble_nk UNIQUE (
        project_id, study_accession, experiment_accession,
        plateid, antigen, family
    )
);

CREATE INDEX IF NOT EXISTS idx_bayes_ensemble_study
    ON madi_results.bayes_ensemble (study_accession, experiment_accession);
CREATE INDEX IF NOT EXISTS idx_bayes_ensemble_curves_fk
    ON madi_results.bayes_ensemble (bayes_curves_id);


-- ─── Done ───────────────────────────────────────────────────────────────────
-- Verify:
--   \dt madi_results.bayes_*
--   SELECT COUNT(*) FROM madi_results.bayes_curves;
