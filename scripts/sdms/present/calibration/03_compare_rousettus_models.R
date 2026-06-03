#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 03_compare_rousettus_models.R ----
# -----------------------------------------------------------------------------|
# Purpose: Compare Gonzalo's saved Rousettus aegyptiacus SDM with regenerated
#          full runs from the present-day calibration workflow.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dismo)
})

# -----------------------------------------------------------------------------|
# Paths and model inputs ----
# -----------------------------------------------------------------------------|

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(
  repo,
  "sdms",
  "runs",
  "chikungunya",
  "calibration",
  "regenerated_models",
  "Rousettus_aegyptiacus",
  "comparison_to_gonzalo"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_paths <- c(
  gonzalo = file.path(repo, "sdms", "models", "Rousettus aegyptiacus", "Rousettus aegyptiacus.rds"),
  our_all_layers = file.path(
    repo,
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "regenerated_models",
    "Rousettus_aegyptiacus",
    "Rousettus_aegyptiacus__spatial-spp__iucn_complete_all__all__2015_2025__bk8000__feature_grid_beta4-8__select10__threads4__model.rds"
  ),
  our_bio_elev = file.path(
    repo,
    "sdms",
    "runs",
    "chikungunya",
    "calibration",
    "regenerated_models",
    "Rousettus_aegyptiacus",
    "Rousettus_aegyptiacus__spatial-spp__iucn_complete_all__bio-elev__2015_2025__bk8000__feature_grid_beta4-8__select10__threads4__model.rds"
  )
)

missing_models <- model_paths[!file.exists(model_paths)]
if (length(missing_models) > 0) {
  stop("Missing model files: ", paste(missing_models, collapse = ", "), call. = FALSE)
}

models <- lapply(model_paths, readRDS)

# -----------------------------------------------------------------------------|
# Fit-performance summaries ----
# -----------------------------------------------------------------------------|

summarise_performance <- function(name, model) {
  params <- model$params
  aicc <- model$AICc
  best <- aicc[which.min(aicc$delta.AICc), , drop = FALSE]
  best_model_col <- if ("model.id" %in% names(best)) "model.id" else "mod"

  data.frame(
    model = name,
    n_retained_models = length(model$mods),
    n_presence = unique(params$n_presence)[1],
    n_background = unique(params$n_background)[1],
    n_variables = length(model$variables),
    variables = paste(model$variables, collapse = "; "),
    train_auc_mean = mean(params$train.AUC.MAXENT, na.rm = TRUE),
    test_auc_mean = mean(params$test.AUC.MAXENT, na.rm = TRUE),
    diff_auc_mean = mean(params$diff.AUC, na.rm = TRUE),
    tss_mean_test_mean = mean(params$TSS.mean.TEST, na.rm = TRUE),
    tss_max_test_mean = mean(params$TSS.max.TEST, na.rm = TRUE),
    boyce_mean = mean(params$C.Boyce, na.rm = TRUE),
    threshold_maxent_median = median(params$TSS.threshold.MAXENT, na.rm = TRUE),
    threshold_test_median = median(params$TSS.threshold.TEST, na.rm = TRUE),
    best_aicc_model = best[[best_model_col]][1],
    best_aicc_ncoefs = best$ncoefs[1],
    best_aicc = best$AICc[1],
    best_aicc_weight = best$w.AIC[1],
    stringsAsFactors = FALSE
  )
}

perf_summary <- do.call(
  rbind,
  Map(summarise_performance, names(models), models)
)
write.csv(
  perf_summary,
  file.path(out_dir, "performance_summary.csv"),
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------|
# Prediction grid over Gonzalo's study area ----
# -----------------------------------------------------------------------------|

union_vars <- unique(unlist(lapply(models, function(model) model$variables)))
predictor_path <- file.path(repo, "sdms", "cache", "Resample_rast.tif")
predictors <- terra::rast(predictor_path)[[union_vars]]

gonzalo_area <- sf::st_as_sf(
  data.frame(name = "gonzalo_study_area", geometry = models$gonzalo$study.area)
)
sf::st_crs(gonzalo_area) <- 4326
gonzalo_area_vect <- terra::vect(gonzalo_area)

predictor_area <- terra::crop(predictors, gonzalo_area_vect)
predictor_area <- terra::mask(predictor_area, gonzalo_area_vect)
predictor_area_5min <- terra::aggregate(predictor_area, fact = 2, fun = mean, na.rm = TRUE)
predictor_area_5min <- terra::mask(predictor_area_5min, gonzalo_area_vect)

values_df <- terra::as.data.frame(predictor_area_5min, cells = TRUE, na.rm = FALSE)
template <- predictor_area_5min[[1]]

predict_average <- function(model, values_df, template) {
  # Recompute predictions from the saved MaxEnt objects because saved raster
  # fields can contain stale terra external pointers after readRDS().
  vars <- model$variables
  complete_idx <- complete.cases(values_df[, vars, drop = FALSE])
  output <- rep(NA_real_, nrow(values_df))

  if (any(complete_idx)) {
    dat <- values_df[complete_idx, vars, drop = FALSE]
    pred_mat <- vapply(
      model$mods,
      function(mod) as.numeric(predict(mod, dat, args = "cloglog")),
      numeric(nrow(dat))
    )
    output[complete_idx] <- rowMeans(pred_mat, na.rm = TRUE)
  }

  raster <- template
  terra::values(raster) <- output
  raster
}

prediction_rasters <- lapply(
  models,
  predict_average,
  values_df = values_df,
  template = template
)
names(prediction_rasters) <- names(models)

pred_stack <- terra::rast(prediction_rasters)
names(pred_stack) <- names(models)
terra::writeRaster(
  pred_stack,
  file.path(out_dir, "average_predictions_gonzalo_extent_5min.tif"),
  overwrite = TRUE
)

prediction_summary <- do.call(
  rbind,
  lapply(names(prediction_rasters), function(name) {
    raster <- prediction_rasters[[name]]
    global_stats <- terra::global(raster, c("min", "max", "mean", "sd"), na.rm = TRUE)
    data.frame(
      model = name,
      prediction_cells = terra::global(!is.na(raster), "sum", na.rm = TRUE)[1, 1],
      min = global_stats[1, "min"],
      max = global_stats[1, "max"],
      mean = global_stats[1, "mean"],
      sd = global_stats[1, "sd"],
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  prediction_summary,
  file.path(out_dir, "prediction_surface_summary_5min.csv"),
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------|
# Continuous and thresholded overlap summaries ----
# -----------------------------------------------------------------------------|

threshold_for <- function(model_name) {
  perf_summary$threshold_maxent_median[perf_summary$model == model_name]
}

cell_area <- terra::cellSize(template, unit = "km")
comparison_rows <- list()
threshold_rows <- list()

for (model_name in setdiff(names(prediction_rasters), "gonzalo")) {
  reference <- prediction_rasters$gonzalo
  comparison <- prediction_rasters[[model_name]]
  values <- terra::values(c(reference, comparison), mat = FALSE, dataframe = TRUE)
  names(values) <- c("gonzalo", model_name)
  complete_idx <- complete.cases(values)
  pair_values <- values[complete_idx, , drop = FALSE]
  difference <- pair_values[[model_name]] - pair_values$gonzalo

  comparison_rows[[model_name]] <- data.frame(
    comparison = paste(model_name, "vs gonzalo"),
    common_cells = nrow(pair_values),
    pearson = cor(pair_values[[model_name]], pair_values$gonzalo, method = "pearson"),
    spearman = cor(pair_values[[model_name]], pair_values$gonzalo, method = "spearman"),
    rmse = sqrt(mean(difference^2)),
    mae = mean(abs(difference)),
    mean_difference = mean(difference),
    median_difference = median(difference),
    q05_difference = as.numeric(quantile(difference, 0.05)),
    q95_difference = as.numeric(quantile(difference, 0.95)),
    stringsAsFactors = FALSE
  )

  reference_suitable <- reference >= threshold_for("gonzalo")
  comparison_suitable <- comparison >= threshold_for(model_name)
  both_suitable <- reference_suitable & comparison_suitable
  either_suitable <- reference_suitable | comparison_suitable
  reference_only <- reference_suitable & !comparison_suitable
  comparison_only <- comparison_suitable & !reference_suitable
  valid <- !is.na(reference_suitable) & !is.na(comparison_suitable)
  valid_values <- terra::values(valid, mat = FALSE) == 1
  area_values <- terra::values(cell_area, mat = FALSE)

  area_for <- function(mask) {
    mask_values <- terra::values(mask, mat = FALSE)
    sum(area_values[mask_values == 1 & valid_values], na.rm = TRUE)
  }

  reference_area <- area_for(reference_suitable)
  comparison_area <- area_for(comparison_suitable)
  both_area <- area_for(both_suitable)
  either_area <- area_for(either_suitable)

  threshold_rows[[model_name]] <- data.frame(
    comparison = paste(model_name, "vs gonzalo"),
    gonzalo_threshold = threshold_for("gonzalo"),
    comparison_threshold = threshold_for(model_name),
    gonzalo_suitable_area_km2 = reference_area,
    comparison_suitable_area_km2 = comparison_area,
    shared_suitable_area_km2 = both_area,
    union_suitable_area_km2 = either_area,
    gonzalo_only_suitable_area_km2 = area_for(reference_only),
    comparison_only_suitable_area_km2 = area_for(comparison_only),
    jaccard_suitable_area = both_area / either_area,
    gonzalo_suitable_captured = both_area / reference_area,
    comparison_suitable_inside_gonzalo = both_area / comparison_area,
    stringsAsFactors = FALSE
  )
}

comparison_summary <- do.call(rbind, comparison_rows)
threshold_summary <- do.call(rbind, threshold_rows)
write.csv(
  comparison_summary,
  file.path(out_dir, "continuous_prediction_comparison_5min.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  threshold_summary,
  file.path(out_dir, "thresholded_prediction_overlap_5min.csv"),
  row.names = FALSE,
  na = ""
)

diff_stack <- terra::rast(list(
  bio_elev_minus_gonzalo = prediction_rasters$our_bio_elev - prediction_rasters$gonzalo,
  all_layers_minus_gonzalo = prediction_rasters$our_all_layers - prediction_rasters$gonzalo
))
terra::writeRaster(
  diff_stack,
  file.path(out_dir, "prediction_differences_minus_gonzalo_5min.tif"),
  overwrite = TRUE
)

overlap_bio <- (prediction_rasters$gonzalo >= threshold_for("gonzalo")) * 1 +
  (prediction_rasters$our_bio_elev >= threshold_for("our_bio_elev")) * 2
overlap_all <- (prediction_rasters$gonzalo >= threshold_for("gonzalo")) * 1 +
  (prediction_rasters$our_all_layers >= threshold_for("our_all_layers")) * 2
names(overlap_bio) <- "bio_elev_overlap_code"
names(overlap_all) <- "all_layers_overlap_code"
terra::writeRaster(
  terra::rast(list(overlap_bio, overlap_all)),
  file.path(out_dir, "thresholded_overlap_codes_5min.tif"),
  overwrite = TRUE
)

# -----------------------------------------------------------------------------|
# Diagnostic map and markdown report ----
# -----------------------------------------------------------------------------|

png(file.path(out_dir, "prediction_comparison_maps_5min.png"), width = 2400, height = 1600, res = 180)
par(mfrow = c(2, 3), mar = c(2, 2, 3, 5))
plot(prediction_rasters$gonzalo, main = "Gonzalo avg prediction", col = hcl.colors(100, "YlGnBu"))
plot(prediction_rasters$our_bio_elev, main = "Our bio-elev avg prediction", col = hcl.colors(100, "YlGnBu"))
plot(
  diff_stack$bio_elev_minus_gonzalo,
  main = "Bio-elev - Gonzalo",
  col = hcl.colors(100, "RdBu", rev = TRUE),
  range = c(-1, 1)
)
plot(prediction_rasters$our_all_layers, main = "Our all-layer avg prediction", col = hcl.colors(100, "YlGnBu"))
plot(
  diff_stack$all_layers_minus_gonzalo,
  main = "All-layer - Gonzalo",
  col = hcl.colors(100, "RdBu", rev = TRUE),
  range = c(-1, 1)
)
plot(
  overlap_bio,
  main = "Threshold overlap: bio-elev vs Gonzalo",
  col = c("grey90", "#d95f02", "#1b9e77", "#4c3b8f"),
  plg = list(title = "0 none, 1 Gonzalo, 2 ours, 3 both")
)
dev.off()

fmt <- function(x, digits = 3) {
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

md_path <- file.path(out_dir, "Rousettus_aegyptiacus_model_comparison.md")
report_lines <- c(
  "# Rousettus aegyptiacus model comparison",
  "",
  "Compared Gonzalo's saved model with two regenerated full runs:",
  "",
  "- `our_bio_elev`: 20 WorldClim bioclim/elevation predictors passed to AutoMaxent.",
  "- `our_all_layers`: full 26-layer stack passed to AutoMaxent.",
  "",
  "Saved `SpatRaster` prediction fields in the `.rds` objects had stale terra external pointers after `readRDS()`, so prediction surfaces were regenerated from the saved MaxEnt model objects and the predictor stack. Surface comparisons use Gonzalo's saved study area and a 5 arc-minute diagnostic grid.",
  "",
  "## Fit and Performance",
  "",
  paste(
    capture.output(print(
      perf_summary[, c(
        "model",
        "n_presence",
        "n_variables",
        "test_auc_mean",
        "tss_mean_test_mean",
        "tss_max_test_mean",
        "boyce_mean",
        "threshold_maxent_median"
      )],
      row.names = FALSE
    )),
    collapse = "\n"
  ),
  "",
  "## Continuous Prediction Comparison",
  "",
  paste(capture.output(print(comparison_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Thresholded Suitable-Area Overlap",
  "",
  paste(capture.output(print(threshold_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Main Read",
  "",
  "- The bio-elev regenerated model is much closer to Gonzalo's predictor set: it keeps all 8 Gonzalo variables and adds only `bio_3`.",
  paste0(
    "- Bio-elev has ",
    perf_summary$n_presence[perf_summary$model == "our_bio_elev"],
    " presences versus Gonzalo's ",
    perf_summary$n_presence[perf_summary$model == "gonzalo"],
    "."
  ),
  paste0(
    "- Bio-elev's mean test AUC is ",
    fmt(perf_summary$test_auc_mean[perf_summary$model == "our_bio_elev"]),
    " versus Gonzalo's ",
    fmt(perf_summary$test_auc_mean[perf_summary$model == "gonzalo"]),
    "."
  ),
  "- All-layer has higher test AUC/TSS in this saved run, but it uses fewer presences and includes predictors that do not appear in Gonzalo's broader saved model set, so it is less useful as a reproduction target.",
  "",
  "## Output Files",
  "",
  "- `performance_summary.csv`",
  "- `prediction_surface_summary_5min.csv`",
  "- `continuous_prediction_comparison_5min.csv`",
  "- `thresholded_prediction_overlap_5min.csv`",
  "- `average_predictions_gonzalo_extent_5min.tif`",
  "- `prediction_differences_minus_gonzalo_5min.tif`",
  "- `thresholded_overlap_codes_5min.tif`",
  "- `prediction_comparison_maps_5min.png`"
)
writeLines(report_lines, md_path)

cat("out_dir=", out_dir, "\n", sep = "")
cat("report=", md_path, "\n", sep = "")
print(
  perf_summary[, c(
    "model",
    "n_presence",
    "n_variables",
    "test_auc_mean",
    "tss_mean_test_mean",
    "tss_max_test_mean",
    "boyce_mean"
  )],
  row.names = FALSE
)
print(comparison_summary, row.names = FALSE)
print(threshold_summary, row.names = FALSE)
