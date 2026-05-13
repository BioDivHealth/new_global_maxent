library(dplyr)
library(readr)
library(tools)

source(here::here("scripts", "associations", "who_don_v2", "who_don_v2_io.R"))

v2_csv_shape <- function(path) {
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x)) {
    return(tibble::tibble(row_count = NA_integer_, column_count = NA_integer_))
  }
  tibble::tibble(row_count = nrow(x), column_count = ncol(x))
}

v2_seeded_baseline_manifest <- function() {
  root <- who_don_v2_output_dir()
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir == FALSE]
  files <- files[!grepl("/qa/v2_seeded_baseline_", files)]

  manifest <- tibble::tibble(
    path = files,
    relative_path = sub(paste0("^", root, "/?"), "", files),
    file_size_bytes = file.info(files)$size,
    md5 = as.character(tools::md5sum(files)),
    file_type = tools::file_ext(files)
  ) %>%
    arrange(relative_path)

  csv_shapes <- manifest %>%
    filter(file_type == "csv") %>%
    rowwise() %>%
    mutate(shape = list(v2_csv_shape(path))) %>%
    tidyr::unnest(shape) %>%
    ungroup() %>%
    select(relative_path, row_count, column_count)

  manifest %>%
    left_join(csv_shapes, by = "relative_path") %>%
    select(relative_path, file_type, file_size_bytes, md5, row_count, column_count)
}

v2_write_seeded_baseline <- function() {
  who_don_v2_ensure_dirs()

  manifest <- v2_seeded_baseline_manifest()
  counts <- manifest %>%
    filter(file_type == "csv") %>%
    transmute(relative_path, row_count, column_count)

  v2_write_csv(manifest, who_don_v2_output_dir("qa", "v2_seeded_baseline_manifest.csv"))
  v2_write_csv(counts, who_don_v2_output_dir("qa", "v2_seeded_baseline_counts.csv"))

  invisible(list(manifest = manifest, counts = counts))
}

