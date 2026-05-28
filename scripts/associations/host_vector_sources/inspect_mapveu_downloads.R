# ------------------------------------------------------------------------------
# Inspect MapVEu text-table downloads
# ------------------------------------------------------------------------------

library(pacman)
p_load(here, readr)

source(here("scripts", "associations", "working_inputs.R"))

mapveu_files <- list.files(
  mapveu_dir,
  pattern = "\\.(txt|tsv|csv)$",
  full.names = TRUE
)

if (length(mapveu_files) == 0) {
  stop("No tabular files found in: ", mapveu_dir)
}

cat("MapVEu directory:", mapveu_dir, "\n")
cat("Files found:", length(mapveu_files), "\n\n")

for (path in sort(mapveu_files)) {
  cat("============================================================\n")
  cat("FILE:", basename(path), "\n")

  preview <- read_tsv(
    path,
    n_max = 5,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA")
  )

  cat("Columns:", ncol(preview), "\n")
  cat("Column names:\n")
  print(names(preview))

  cat("\nHead:\n")
  print(preview)
  cat("\n")
}
