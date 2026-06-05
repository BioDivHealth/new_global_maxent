#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------|
# 07_plot_combined_vector_occurrence_maps.R ----
# -----------------------------------------------------------------------------|
# Purpose: Plot cleaned combined Chikungunya vector occurrences by source.
#
# Maps are diagnostics for the model-facing occurrence layer. They read cleaned
# combined files and do not modify GBIF, VectorMap, MapVEu, or combined inputs.
# -----------------------------------------------------------------------------|

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required.", call. = FALSE)
  }
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package `maps` is required.", call. = FALSE)
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package `grid` is required.", call. = FALSE)
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package `here` is required.", call. = FALSE)
  }
})

source(file.path(here::here(), "scripts", "sdms", "present", "utils.R"))

# -----------------------------------------------------------------------------|
# RStudio config: edit this block before sourcing the script ----
# -----------------------------------------------------------------------------|

if (!exists("batch_config", inherits = FALSE)) {
  batch_config <- list(
    roles = "vector",
    occurrence_method = "combined",
    start_year = 1970,
    end_year = as.integer(format(Sys.Date(), "%Y")),
    include_missing_year_records = TRUE,
    species_filter = character(),
    max_species = Inf,
    point_size = 0.75,
    point_alpha = 0.65,
    adaptive_point_size = TRUE,
    zoom_buffer_fraction = 0.08,
    zoom_min_buffer_degrees = 2,
    zoom_min_width_degrees = 10,
    zoom_min_height_degrees = 8,
    timestamped_run_dir = FALSE,
    run_label = NA_character_,
    dry_run = FALSE
  )
}

# -----------------------------------------------------------------------------|
# Internal defaults ----
# -----------------------------------------------------------------------------|

default_batch_config <- list(
  target_manifest_path = file.path(repo_root(), "sdms", "runs", "chikungunya", "sdm_target_manifest.csv"),
  occurrence_root = file.path(repo_root(), "sdms", "runs", "vector_sdm_push", "occurrences"),
  map_run_root = file.path(repo_root(), "sdms", "runs", "chikungunya", "maps", "combined_vector_occurrences"),
  occurrence_method = "combined",
  roles = "vector",
  include_not_needed = FALSE,
  include_already_available = FALSE,
  species_filter = character(),
  max_species = Inf,
  start_year = 1970,
  end_year = as.integer(format(Sys.Date(), "%Y")),
  include_missing_year_records = TRUE,
  point_size = 0.75,
  point_alpha = 0.65,
  adaptive_point_size = TRUE,
  sparse_point_size = 1.8,
  medium_point_size = 1.25,
  dense_point_size = 0.95,
  zoom_buffer_fraction = 0.08,
  zoom_min_buffer_degrees = 2,
  zoom_min_width_degrees = 10,
  zoom_min_height_degrees = 8,
  map_width = 12,
  map_height = 6.5,
  map_dpi = 180,
  timestamped_run_dir = FALSE,
  run_label = NA_character_,
  dry_run = FALSE
)

batch_config <- utils::modifyList(default_batch_config, batch_config)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# -----------------------------------------------------------------------------|
# Config helpers ----
# -----------------------------------------------------------------------------|

config_arg <- function(key, config_key = gsub("-", "_", key)) {
  get_arg(args, key, batch_config[[config_key]])
}

target_manifest_path <- config_arg("target-manifest-path")
occurrence_root <- config_arg("occurrence-root")
map_run_root <- config_arg("map-run-root")
occurrence_method <- config_arg("occurrence-method")
roles <- split_arg(config_arg("roles"))
include_not_needed <- as_logical_arg(config_arg("include-not-needed"))
include_already_available <- as_logical_arg(config_arg("include-already-available"))
species_filter <- split_arg(config_arg("species-filter"))
max_species <- as.numeric(config_arg("max-species"))
start_year <- as.integer(config_arg("start-year"))
end_year <- as.integer(config_arg("end-year"))
include_missing_year_records <- as_logical_arg(config_arg("include-missing-year-records"))
point_size <- as.numeric(config_arg("point-size"))
point_alpha <- as.numeric(config_arg("point-alpha"))
adaptive_point_size <- as_logical_arg(config_arg("adaptive-point-size"))
sparse_point_size <- as.numeric(config_arg("sparse-point-size"))
medium_point_size <- as.numeric(config_arg("medium-point-size"))
dense_point_size <- as.numeric(config_arg("dense-point-size"))
zoom_buffer_fraction <- as.numeric(config_arg("zoom-buffer-fraction"))
zoom_min_buffer_degrees <- as.numeric(config_arg("zoom-min-buffer-degrees"))
zoom_min_width_degrees <- as.numeric(config_arg("zoom-min-width-degrees"))
zoom_min_height_degrees <- as.numeric(config_arg("zoom-min-height-degrees"))
map_width <- as.numeric(config_arg("map-width"))
map_height <- as.numeric(config_arg("map-height"))
map_dpi <- as.integer(config_arg("map-dpi"))
timestamped_run_dir <- as_logical_arg(config_arg("timestamped-run-dir"))
run_label <- config_arg("run-label")
dry_run <- as_logical_arg(config_arg("dry-run")) || has_flag(args, "dry-run")

if (!file.exists(target_manifest_path)) {
  stop("Missing Chikungunya SDM target manifest: ", target_manifest_path, call. = FALSE)
}

# -----------------------------------------------------------------------------|
# Helpers ----
# -----------------------------------------------------------------------------|

utc_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

default_run_label <- function() {
  missing_year_suffix <- if (include_missing_year_records) "_including_missing_years" else "_known_years_only"
  paste0(occurrence_method, "_", start_year, "_", end_year, missing_year_suffix)
}

cleaned_occurrence_path <- function(species) {
  species_safe <- safe_species_name(species)
  file.path(
    occurrence_root,
    species_safe,
    occurrence_method,
    "cleaned",
    paste0(species_safe, "_cleaned.csv")
  )
}

source_plot_class <- function(data) {
  has_gbif <- data$gbif_download_row_count > 0
  has_vectormap <- data$vectormap_row_count > 0
  has_mapveu <- data$mapveu_row_count > 0

  out <- rep(NA_character_, nrow(data))
  out[has_gbif & !has_vectormap & !has_mapveu] <- "GBIF only"
  out[!has_gbif & has_vectormap & !has_mapveu] <- "VectorMap only"
  out[!has_gbif & !has_vectormap & has_mapveu] <- "MapVEu only"
  out[has_gbif & has_vectormap & !has_mapveu] <- "GBIF + VectorMap"
  out[has_gbif & !has_vectormap & has_mapveu] <- "GBIF + MapVEu"
  out[!has_gbif & has_vectormap & has_mapveu] <- "VectorMap + MapVEu"
  out[has_gbif & has_vectormap & has_mapveu] <- "GBIF + VectorMap + MapVEu"
  out
}

read_cleaned_occurrences <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required_cols <- c("decimalLongitude", "decimalLatitude")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Cleaned occurrence file is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  for (col in c("gbif_download_row_count", "vectormap_row_count", "mapveu_row_count")) {
    if (!col %in% names(data)) {
      data[[col]] <- 0L
    }
    data[[col]] <- suppressWarnings(as.integer(data[[col]]))
    data[[col]][is.na(data[[col]])] <- 0L
  }

  data$source_plot_class <- source_plot_class(data)
  data
}

filter_map_year_window <- function(data) {
  for (col in c("year_min", "year_max", "year")) {
    if (!col %in% names(data)) {
      data[[col]] <- NA_integer_
    }
    data[[col]] <- suppressWarnings(as.integer(data[[col]]))
  }

  year_min <- ifelse(is.na(data$year_min), data$year, data$year_min)
  year_max <- ifelse(is.na(data$year_max), data$year, data$year_max)
  has_known_year <- !is.na(year_min) | !is.na(year_max)
  overlaps_window <- has_known_year &
    (is.na(year_max) | year_max >= start_year) &
    (is.na(year_min) | year_min <= end_year)

  keep <- overlaps_window | (!has_known_year & include_missing_year_records)
  data[keep, , drop = FALSE]
}

plot_point_size <- function(row_count) {
  if (!adaptive_point_size) {
    return(point_size)
  }

  if (row_count <= 75) {
    return(sparse_point_size)
  }
  if (row_count <= 750) {
    return(medium_point_size)
  }
  if (row_count <= 5000) {
    return(dense_point_size)
  }

  point_size
}

clamp <- function(x, lower, upper) {
  min(max(x, lower), upper)
}

expand_bounds <- function(min_value, max_value, lower_limit, upper_limit, min_span) {
  span <- max_value - min_value
  if (is.na(span) || span < min_span) {
    center <- mean(c(min_value, max_value), na.rm = TRUE)
    min_value <- center - min_span / 2
    max_value <- center + min_span / 2
  }

  if (min_value < lower_limit) {
    max_value <- max_value + (lower_limit - min_value)
    min_value <- lower_limit
  }
  if (max_value > upper_limit) {
    min_value <- min_value - (max_value - upper_limit)
    max_value <- upper_limit
  }

  c(
    min = clamp(min_value, lower_limit, upper_limit),
    max = clamp(max_value, lower_limit, upper_limit)
  )
}

occurrence_bounds <- function(data) {
  lon <- data$decimalLongitude[is.finite(data$decimalLongitude)]
  lat <- data$decimalLatitude[is.finite(data$decimalLatitude)]
  if (length(lon) == 0 || length(lat) == 0) {
    return(c(xmin = -180, xmax = 180, ymin = -60, ymax = 85))
  }

  lon_span <- diff(range(lon))
  lat_span <- diff(range(lat))
  lon_buffer <- max(lon_span * zoom_buffer_fraction, zoom_min_buffer_degrees)
  lat_buffer <- max(lat_span * zoom_buffer_fraction, zoom_min_buffer_degrees)

  x_bounds <- expand_bounds(
    min(lon) - lon_buffer,
    max(lon) + lon_buffer,
    -180,
    180,
    zoom_min_width_degrees
  )
  y_bounds <- expand_bounds(
    min(lat) - lat_buffer,
    max(lat) + lat_buffer,
    -60,
    85,
    zoom_min_height_degrees
  )

  c(xmin = x_bounds[["min"]], xmax = x_bounds[["max"]], ymin = y_bounds[["min"]], ymax = y_bounds[["max"]])
}

extract_legend <- function(plot) {
  grob <- ggplot2::ggplotGrob(plot)
  guide_index <- which(vapply(grob$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(guide_index) == 0) {
    return(NULL)
  }

  grob$grobs[[guide_index[[1]]]]
}

save_map_layout <- function(main_plot, inset_plot, legend_grob, output_path) {
  grDevices::png(
    filename = output_path,
    width = map_width * map_dpi,
    height = map_height * map_dpi,
    res = map_dpi
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2,
    ncol = 2,
    widths = grid::unit(c(4.9, 1.25), "null"),
    heights = grid::unit(c(1, 0.18), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  print(
    main_plot + ggplot2::theme(legend.position = "none"),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  print(
    inset_plot,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2)
  )

  if (!is.null(legend_grob)) {
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
    grid::grid.draw(legend_grob)
    grid::popViewport()
  }

  grid::popViewport()
}

world_overview_plot <- function(world, bounds) {
  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "#F2EFE9",
      color = "#C8C4BA",
      linewidth = 0.1
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = bounds[["xmin"]],
        xmax = bounds[["xmax"]],
        ymin = bounds[["ymin"]],
        ymax = bounds[["ymax"]]
      ),
      inherit.aes = FALSE,
      fill = NA,
      color = "#1B1B1B",
      linewidth = 0.6
    ) +
    ggplot2::coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
    ggplot2::labs(title = "Global context") +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = 0),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = "#B8B8B8", linewidth = 0.35),
      plot.margin = ggplot2::margin(28, 8, 8, 8)
    )
}

plot_species_map <- function(data, species, output_path, bounds) {
  world <- ggplot2::map_data("world")
  effective_point_size <- plot_point_size(nrow(data))
  plot_order <- c(
    "GBIF only",
    "VectorMap only",
    "MapVEu only",
    "GBIF + VectorMap",
    "GBIF + MapVEu",
    "VectorMap + MapVEu",
    "GBIF + VectorMap + MapVEu"
  )
  data$source_plot_class <- factor(data$source_plot_class, levels = plot_order)
  legend_breaks <- plot_order[plot_order %in% as.character(unique(data$source_plot_class))]

  palette <- c(
    "GBIF only" = "#2B6CB0",
    "VectorMap only" = "#C43C39",
    "MapVEu only" = "#C99A13",
    "GBIF + VectorMap" = "#6B3FA0",
    "GBIF + MapVEu" = "#2F8F57",
    "VectorMap + MapVEu" = "#E07A1F",
    "GBIF + VectorMap + MapVEu" = "#3B2F2F"
  )

  title <- paste0(species, " cleaned combined occurrences")
  subtitle <- paste0(
    start_year,
    "-",
    end_year,
    "; n = ",
    format(nrow(data), big.mark = ","),
    " cleaned coordinate rows"
  )

  main_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = world,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "#F2EFE9",
      color = "#C8C4BA",
      linewidth = 0.15
    ) +
    ggplot2::geom_point(
      data = data,
      ggplot2::aes(
        x = decimalLongitude,
        y = decimalLatitude,
        color = source_plot_class
      ),
      size = effective_point_size,
      alpha = point_alpha,
      stroke = 0
    ) +
    ggplot2::scale_color_manual(
      values = palette,
      breaks = legend_breaks,
      drop = TRUE,
      na.translate = FALSE,
      name = "Source"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4.2, alpha = 1),
        nrow = 1,
        byrow = TRUE
      )
    ) +
    ggplot2::coord_quickmap(
      xlim = c(bounds[["xmin"]], bounds[["xmax"]]),
      ylim = c(bounds[["ymin"]], bounds[["ymax"]]),
      expand = FALSE
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_line(color = "#E4E1DA", linewidth = 0.25),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "#555555"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width = grid::unit(1.2, "lines"),
      legend.spacing.x = grid::unit(0.8, "lines"),
      axis.text = ggplot2::element_text(color = "#555555")
    )

  legend_grob <- extract_legend(main_plot)
  save_map_layout(main_plot, world_overview_plot(world, bounds), legend_grob, output_path)
}

# -----------------------------------------------------------------------------|
# Select targets and plot maps ----
# -----------------------------------------------------------------------------|

target_manifest <- read.csv(target_manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
targets <- select_sdm_targets(
  target_manifest = target_manifest,
  roles = roles,
  species_filter = species_filter,
  include_not_needed = include_not_needed,
  include_already_available = include_already_available,
  max_species = max_species
)

timestamp <- paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "_pid", Sys.getpid())
if (is.na(run_label) || !nzchar(run_label)) {
  run_label <- default_run_label()
}
run_dir <- file.path(map_run_root, if (timestamped_run_dir) timestamp else run_label)
map_dir <- file.path(run_dir, "maps")
summary_path <- file.path(run_dir, "combined_vector_occurrence_map_summary.csv")

cat("Selected target species:", nrow(targets), "\n")
if (dry_run) {
  cat("Dry run: maps and summaries will not be written.\n")
}

rows <- vector("list", nrow(targets))
for (i in seq_len(nrow(targets))) {
  species <- targets$species_name_canonical[[i]]
  species_safe <- safe_species_name(species)
  cleaned_path <- cleaned_occurrence_path(species)
  map_path <- file.path(map_dir, paste0(species_safe, "_combined_cleaned_sources_map.png"))

  result <- tryCatch(
    {
      data <- read_cleaned_occurrences(cleaned_path)
      if (is.null(data) || nrow(data) == 0) {
        result_row <- data.frame(
          species_name = species,
          map_status = "no_cleaned_records",
          cleaned_path = cleaned_path,
          map_path = NA_character_,
          cleaned_rows = 0,
          mapped_rows = 0,
          year_filtered_rows = 0,
          point_size_used = NA_real_,
          gbif_only_rows = 0,
          vectormap_only_rows = 0,
          mapveu_only_rows = 0,
          gbif_vectormap_rows = 0,
          gbif_mapveu_rows = 0,
          vectormap_mapveu_rows = 0,
          all_sources_rows = 0,
          map_xlim_min = NA_real_,
          map_xlim_max = NA_real_,
          map_ylim_min = NA_real_,
          map_ylim_max = NA_real_,
          prepared_at = utc_now(),
          stringsAsFactors = FALSE
        )
      } else {
        source_cleaned_rows <- nrow(data)
        data <- filter_map_year_window(data)
        year_filtered_rows <- source_cleaned_rows - nrow(data)

        if (nrow(data) == 0) {
          result_row <- data.frame(
            species_name = species,
            map_status = "no_records_after_year_filter",
            cleaned_path = cleaned_path,
            map_path = NA_character_,
            cleaned_rows = source_cleaned_rows,
            mapped_rows = 0,
            year_filtered_rows = year_filtered_rows,
            point_size_used = NA_real_,
            gbif_only_rows = 0,
            vectormap_only_rows = 0,
            mapveu_only_rows = 0,
            gbif_vectormap_rows = 0,
            gbif_mapveu_rows = 0,
            vectormap_mapveu_rows = 0,
            all_sources_rows = 0,
            map_xlim_min = NA_real_,
            map_xlim_max = NA_real_,
            map_ylim_min = NA_real_,
            map_ylim_max = NA_real_,
            prepared_at = utc_now(),
            stringsAsFactors = FALSE
          )
        } else {
          bounds <- occurrence_bounds(data)
          effective_point_size <- plot_point_size(nrow(data))
          counts <- table(factor(
            data$source_plot_class,
            levels = c(
              "GBIF only",
              "VectorMap only",
              "MapVEu only",
              "GBIF + VectorMap",
              "GBIF + MapVEu",
              "VectorMap + MapVEu",
              "GBIF + VectorMap + MapVEu"
            )
          ))

          if (!dry_run) {
            ensure_dir(map_dir)
            plot_species_map(data, species, map_path, bounds)
          }

          result_row <- data.frame(
            species_name = species,
            map_status = if (!dry_run) "mapped" else "dry_run_ready",
            cleaned_path = cleaned_path,
            map_path = if (!dry_run && file.exists(map_path)) map_path else NA_character_,
            cleaned_rows = source_cleaned_rows,
            mapped_rows = nrow(data),
            year_filtered_rows = year_filtered_rows,
            point_size_used = effective_point_size,
            gbif_only_rows = as.integer(counts[["GBIF only"]]),
            vectormap_only_rows = as.integer(counts[["VectorMap only"]]),
            mapveu_only_rows = as.integer(counts[["MapVEu only"]]),
            gbif_vectormap_rows = as.integer(counts[["GBIF + VectorMap"]]),
            gbif_mapveu_rows = as.integer(counts[["GBIF + MapVEu"]]),
            vectormap_mapveu_rows = as.integer(counts[["VectorMap + MapVEu"]]),
            all_sources_rows = as.integer(counts[["GBIF + VectorMap + MapVEu"]]),
            map_xlim_min = bounds[["xmin"]],
            map_xlim_max = bounds[["xmax"]],
            map_ylim_min = bounds[["ymin"]],
            map_ylim_max = bounds[["ymax"]],
            prepared_at = utc_now(),
            stringsAsFactors = FALSE
          )
        }
      }
      result_row
    },
    error = function(err) {
      data.frame(
        species_name = species,
        map_status = "failed_error",
        cleaned_path = cleaned_path,
        map_path = NA_character_,
        cleaned_rows = NA_integer_,
        mapped_rows = NA_integer_,
        year_filtered_rows = NA_integer_,
        point_size_used = NA_real_,
        gbif_only_rows = NA_integer_,
        vectormap_only_rows = NA_integer_,
        mapveu_only_rows = NA_integer_,
        gbif_vectormap_rows = NA_integer_,
        gbif_mapveu_rows = NA_integer_,
        vectormap_mapveu_rows = NA_integer_,
        all_sources_rows = NA_integer_,
        map_xlim_min = NA_real_,
        map_xlim_max = NA_real_,
        map_ylim_min = NA_real_,
        map_ylim_max = NA_real_,
        prepared_at = utc_now(),
        notes = conditionMessage(err),
        stringsAsFactors = FALSE
      )
    }
  )

  rows[[i]] <- result
  if (!dry_run) {
    ensure_dir(run_dir)
    data.table::fwrite(data.table::rbindlist(rows[seq_len(i)], use.names = TRUE, fill = TRUE), summary_path, na = "")
  }

  cat("[", i, "/", nrow(targets), "] ", species, ": ", result$map_status[[1]], "\n", sep = "")
}

summary <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
if (!dry_run) {
  ensure_dir(run_dir)
  data.table::fwrite(summary, summary_path, na = "")
  cat("Wrote map run summary:", summary_path, "\n")
}
