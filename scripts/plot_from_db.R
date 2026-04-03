# =============================================================================
# plot_from_db.R — Reconstruct Bayesian standard curve plot from DB data
#
# Uses pre-computed curve grid (CI ribbon) and CDAN grid (precision profile)
# stored by the worker. Standards are median-aggregated + prozone-corrected
# to match what the fit actually saw.
#
# Usage:  Rscript scripts/plot_from_db.R
# Output: plots/bayes_plate_from_db.html
# =============================================================================

suppressPackageStartupMessages({
  library(RPostgres)
  library(DBI)
  library(dplyr)
  library(plotly)
})

PLATE_ID <- "GAPS_ADCD_PLATE12_06062024"
ANTIGEN  <- "tt"
SOURCE   <- "NIBSC06_140"

conn <- dbConnect(Postgres(),
  dbname = "local_madi_ispi", host = "localhost", port = 5432L,
  user = Sys.getenv("DB_USER", "hardik"), password = Sys.getenv("DB_PASSWORD", ""))

# ── 1. Curve params for this specific source ──────────────────────────────
curve <- dbGetQuery(conn, sprintf(
  "SELECT * FROM madi_results.bayes_curves
   WHERE plateid = '%s' AND antigen = '%s' AND source = '%s'",
  PLATE_ID, ANTIGEN, SOURCE))
# Also look up standard_curve_concentration for the concentration formula
sc_conc_row <- dbGetQuery(conn, sprintf(
  "SELECT DISTINCT standard_curve_concentration FROM madi_results.xmap_antigen_family
   WHERE study_accession = '%s' AND experiment_accession = '%s' AND antigen = '%s'
     AND standard_curve_concentration IS NOT NULL LIMIT 1",
  "MADI_P3_GAPS", "ADCD", ANTIGEN))
sc_conc <- if (nrow(sc_conc_row) > 0) as.numeric(sc_conc_row$standard_curve_concentration[1]) else NULL

stopifnot(nrow(curve) == 1)
nom <- if (!is.null(sc_conc)) sc_conc else as.numeric(curve$nominal_sample_dilution)
cat(sprintf("Using base concentration: %g\n", nom))
cat(sprintf("Curve: %s (global_best=%s)  a=%.1f b=%.2f c=%.4f d=%.1f g=%.3f\n",
            curve$curve_family, curve$global_best_family,
            curve$a, curve$b, curve$c, curve$d, curve$g))

# ── 2. Standards — single source, median-aggregated, prozone-corrected ────
# Each source is fitted independently, so we only show that source's standards.
stds_raw <- dbGetQuery(conn, sprintf(
  "SELECT dilution as dilution_factor, antibody_mfi as mfi
   FROM madi_results.xmap_standard
   WHERE plateid = '%s' AND antigen = '%s' AND source = '%s'
     AND antibody_mfi > 0 AND dilution > 0",
  PLATE_ID, ANTIGEN, SOURCE))

stds_raw$concentration <- nom / stds_raw$dilution_factor

# Median-aggregate per concentration (combines both sources)
stds <- stds_raw |>
  group_by(concentration) |>
  summarise(mfi = median(mfi, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(mfi))

# Apply prozone correction (same as worker_batch.R correct_prozone())
correct_prozone <- function(df, prop_diff = 0.1, dil_scale = 2) {
  max_mfi <- max(df$mfi, na.rm = TRUE)
  conc_at_max <- max(df$concentration[df$mfi == max_mfi])
  post_peak <- df$concentration > conc_at_max
  if (any(post_peak)) {
    df$mfi[post_peak] <- max_mfi +
      (max_mfi - df$mfi[post_peak]) * prop_diff /
      ((df$concentration[post_peak] - conc_at_max) * dil_scale)
  }
  df
}
stds <- correct_prozone(stds)
cat(sprintf("Standards: %d points (prozone-corrected, median-aggregated)\n", nrow(stds)))

# ── 3. Sample predictions for this source ──────────────────────────────────
samples <- dbGetQuery(conn, sprintf(
  "SELECT sampleid, mfi, raw_predicted_concentration as pred_conc,
          conc_lower, conc_upper, pcov, gate_class
   FROM madi_results.bayes_samples
   WHERE plateid = '%s' AND antigen = '%s' AND source = '%s'",
  PLATE_ID, ANTIGEN, SOURCE))
cat(sprintf("Samples: %d\n", nrow(samples)))

# ── 4. Curve grid (CI ribbon from posterior draws) ─────────────────────────
curve_grid <- dbGetQuery(conn, sprintf(
  "SELECT log10_conc, mfi_median, mfi_lower_95, mfi_upper_95
   FROM madi_results.bayes_curve_grid
   WHERE plateid = '%s' AND antigen = '%s' AND source = '%s'
   ORDER BY log10_conc", PLATE_ID, ANTIGEN, SOURCE))
cat(sprintf("Curve grid: %d points\n", nrow(curve_grid)))

# ── 5. CDAN grid (precision profile) ──────────────────────────────────────
cdan_grid <- dbGetQuery(conn, sprintf(
  "SELECT log10_conc, smoothed_cv
   FROM madi_results.bayes_cdan_grid
   WHERE plateid = '%s' AND antigen = '%s' AND source = '%s'
     AND smoothed_cv IS NOT NULL AND smoothed_cv < 60
   ORDER BY log10_conc", PLATE_ID, ANTIGEN, SOURCE))
cat(sprintf("CDAN grid: %d points\n", nrow(cdan_grid)))

dbDisconnect(conn)

# ── 6. Prepare fit data ───────────────────────────────────────────────────
fit_df <- data.frame(
  x       = curve_grid$log10_conc,
  y       = log10(pmax(curve_grid$mfi_median, 1e-9)),
  y_lower = log10(pmax(curve_grid$mfi_lower_95, 1e-9)),
  y_upper = log10(pmax(curve_grid$mfi_upper_95, 1e-9))
)

# ── 7. Colors ─────────────────────────────────────────────────────────────
COL_STD <- "#000000"; COL_FIT <- "#0072B2"; COL_CI <- "rgba(86,180,233,0.20)"
COL_ASYM <- "#999999"; COL_LOQ <- "#D55E00"; COL_RDL <- "#604e97"
COL_SAMP <- "#E69F00"; COL_INFL <- "#009E73"; COL_CDAN <- "#1565C0"

fam <- switch(curve$curve_family, "4pl"="4PL", "5pl"="5PL", "gompertz"="Gompertz",
              curve$curve_family)

y_all <- c(fit_df$y_lower, fit_df$y_upper, log10(pmax(stds$mfi, 1)))
y_lo <- min(y_all, na.rm = TRUE); y_hi <- max(y_all, na.rm = TRUE)

# ── 8. Build plot ──────────────────────────────────────────────────────────
p <- plot_ly() |>
  add_markers(data = stds, x = ~log10(concentration), y = ~log10(mfi),
    name = "Standards", marker = list(color = COL_STD, size = 7),
    hovertemplate = "Conc: %{customdata:.4f}<br>MFI: %{y:.3f}<extra></extra>",
    customdata = stds$concentration) |>
  add_ribbons(data = fit_df, x = ~x, ymin = ~y_lower, ymax = ~y_upper,
    name = "95% CI", fillcolor = COL_CI, line = list(color = "transparent")) |>
  add_lines(data = fit_df, x = ~x, y = ~y,
    name = paste0(fam, " Fit"), line = list(color = COL_FIT, width = 2.5))

# Samples
vs <- samples |> filter(!is.na(pred_conc), pred_conc > 0, !is.na(mfi), mfi > 0)
if (nrow(vs) > 0) {
  vs$hover <- paste0(
    "<b>ID:</b> ", vs$sampleid,
    "<br><b>MFI:</b> ", round(vs$mfi, 1),
    "<br><b>Conc:</b> ", signif(vs$pred_conc, 4),
    "<br><b>95%CI:</b> [", signif(vs$conc_lower, 3), ", ", signif(vs$conc_upper, 3), "]",
    "<br><b>pCoV:</b> ", ifelse(!is.na(vs$pcov), paste0(round(vs$pcov, 1), "%"), "N/A"),
    "<br><b>Gate:</b> ", vs$gate_class)
  p <- p |> add_markers(data = vs, x = ~log10(pred_conc), y = ~log10(mfi),
    name = "Samples", text = ~hover, hoverinfo = "text",
    marker = list(color = COL_SAMP, size = 7, symbol = "diamond-open",
                  line = list(width = 2, color = COL_SAMP)))
}

# Vertical lines helper
add_vline <- function(p, val, nm, lg, col, dash, show) {
  if (!is.na(val) && is.finite(val) && val > 0)
    p |> add_segments(x = log10(val), xend = log10(val), y = y_lo, yend = y_hi,
      name = nm, legendgroup = lg, line = list(color = col, dash = dash, width = 2),
      showlegend = show, hoverinfo = "skip")
  else p
}
p <- p |>
  add_vline(curve$lloq, "LLOQ / ULOQ", "loq", COL_LOQ, "dash", TRUE) |>
  add_vline(curve$uloq, "LLOQ / ULOQ", "loq", COL_LOQ, "dash", FALSE) |>
  add_vline(curve$lrdl, "LRDL / URDL", "rdl", COL_RDL, "dashdot", TRUE) |>
  add_vline(curve$urdl, "LRDL / URDL", "rdl", COL_RDL, "dashdot", FALSE)

# Inflection
if (!is.na(curve$inflect_x) && curve$inflect_x > 0) {
  ix <- log10(curve$inflect_x)
  iy <- approx(fit_df$x, fit_df$y, xout = ix)$y
  p <- p |> add_markers(x = ix, y = iy,
    name = "Inflection (max |dy/dx|)",
    marker = list(color = COL_INFL, size = 10, symbol = "diamond",
                  line = list(width = 1.5, color = "#000000")))
}

# Asymptotes (hidden by default)
if (!is.na(curve$a) && curve$a > 0)
  p <- p |> add_segments(x = min(fit_df$x), xend = max(fit_df$x),
    y = log10(curve$a), yend = log10(curve$a),
    name = "Lower Asymptote", legendgroup = "asym",
    line = list(color = COL_ASYM, width = 1, dash = "dot"), visible = "legendonly")
if (!is.na(curve$d) && curve$d > 0)
  p <- p |> add_segments(x = min(fit_df$x), xend = max(fit_df$x),
    y = log10(curve$d), yend = log10(curve$d),
    name = "Upper Asymptote", legendgroup = "asym",
    line = list(color = COL_ASYM, width = 1, dash = "dot"),
    showlegend = FALSE, visible = "legendonly")

# ── 9. CDAN precision profile (secondary Y axis) ──────────────────────────
if (nrow(cdan_grid) > 0) {
  p <- p |>
    add_lines(data = cdan_grid, x = ~log10_conc, y = ~smoothed_cv,
      name = "Bayesian CDAN Precision Profile", yaxis = "y2",
      line = list(color = COL_CDAN, width = 2.5),
      hovertemplate = "Log10 Conc: %{x:.2f}<br>CV%%: %{y:.1f}%%<extra>CDAN</extra>") |>
    add_lines(x = range(cdan_grid$log10_conc), y = c(20, 20),
      name = "pCoV Threshold: 20%", yaxis = "y2",
      line = list(color = "#e68fac", dash = "dash", width = 1.5), hoverinfo = "skip") |>
    add_lines(x = range(cdan_grid$log10_conc), y = c(15, 15),
      name = "pCoV Threshold: 15%", yaxis = "y2",
      line = list(color = "#4CAF50", dash = "dash", width = 1.5), hoverinfo = "skip")
}

# ── 10. Layout ─────────────────────────────────────────────────────────────
layout_args <- list(
  title = list(text = paste0(fam, " Fit (from DB) \u2014 ", PLATE_ID), font = list(size = 14)),
  xaxis = list(title = "Log\u2081\u2080 Concentration",
    gridcolor = "#E5E5E5", showline = TRUE, linecolor = "#CCCCCC"),
  yaxis = list(title = "Log\u2081\u2080 MFI",
    gridcolor = "#E5E5E5", showline = TRUE, linecolor = "#CCCCCC"),
  plot_bgcolor = "white", paper_bgcolor = "white",
  hovermode = "closest",
  legend = list(x = 1.05, y = 1, xanchor = "left",
    bgcolor = "rgba(255,255,255,0.85)", bordercolor = "#CCCCCC", borderwidth = 1))
if (nrow(cdan_grid) > 0) {
  layout_args$yaxis2 <- list(
    overlaying = "y", side = "right",
    title = "Concentration Uncertainty (pCoV %)",
    range = c(0, 55), showgrid = FALSE, zeroline = FALSE,
    tickfont = list(color = COL_CDAN), titlefont = list(color = COL_CDAN))
}
p <- do.call(plotly::layout, c(list(p), layout_args))

# ── 11. Save ───────────────────────────────────────────────────────────────
dir.create("plots", showWarnings = FALSE)
htmlwidgets::saveWidget(p, "plots/bayes_plate_from_db.html", selfcontained = TRUE)
cat("\nSaved: plots/bayes_plate_from_db.html\n")
