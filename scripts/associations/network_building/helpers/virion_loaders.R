# -----------------------------------------------------------------------------|
# virion_loaders.R ----
# -----------------------------------------------------------------------------|
# Purpose: Reusable VIRION table loader functions without automatically loading
#          the full local VIRION dataset.
# -----------------------------------------------------------------------------|

if (!exists("virion_source_version_dir") || !exists("virion_source_dir")) {
  source(file.path("scripts", "associations", "working_inputs.R"))
}

read_virion_csv <- function(path) {
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) {
      if (grepl("input string .* is invalid", conditionMessage(e), ignore.case = TRUE)) {
        warning(
          "Encountered invalid text encoding while reading ", basename(path),
          "; retrying with Latin1 decoding."
        )
        return(
          readr::read_csv(
            path,
            show_col_types = FALSE,
            locale = readr::locale(encoding = "Latin1")
          )
        )
      }
      stop(e)
    }
  )
}

load_virion_data <- function(
  data_path = virion_source_version_dir,
  files = c(
    "virion.csv.gz",
    "edgelist.csv",
    "taxonomy_host.csv",
    "taxonomy_virus.csv",
    "provenance.csv.gz",
    "detection.csv.gz",
    "temporal.csv.gz"
  )
) {
  if (!dir.exists(data_path)) {
    stop(
      "VIRION data directory not found. Please download data from Zenodo or ",
      "use virionData package."
    )
  }

  data_list <- list()

  for (file in files) {
    file_path <- file.path(data_path, file)

    if (file.exists(file_path)) {
      cat("Loading:", file, "\n")

      table_name <- gsub("\\.csv(\\.gz)?$", "", file)
      data_list[[table_name]] <- read_virion_csv(file_path)
    } else {
      warning("File not found:", file_path)
    }
  }

  data_list
}

load_virion_package <- function(
  version = "latest",
  tables = c(
    "virion",
    "edgelist",
    "taxonomy_host",
    "taxonomy_virus",
    "provenance",
    "detection",
    "temporal"
  )
) {
  if (!require(virionData, quietly = TRUE)) {
    cat("virionData package not installed. Installing...\n")
    if (!require(remotes, quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("viralemergence/virionData")
    library(virionData)
  }

  virion_data <- list()

  cat("Loading VIRION data tables...\n")
  cat("  - Downloading VIRION data from Zenodo...\n")
  data_path <- virionData::get_versioned_data(version = version, dir_path = virion_source_dir)

  if ("virion" %in% tables) {
    cat("  - Loading main virion interactions...\n")
    virion_data$virion <- read_virion_csv(file.path(data_path, "virion.csv.gz"))
  }

  if ("edgelist" %in% tables) {
    cat("  - Loading edgelist...\n")
    virion_data$edgelist <- readr::read_csv(file.path(data_path, "edgelist.csv"), show_col_types = FALSE)
  }

  if ("taxonomy_host" %in% tables) {
    cat("  - Loading host taxonomy...\n")
    virion_data$taxonomy_host <- readr::read_csv(file.path(data_path, "taxonomy_host.csv"), show_col_types = FALSE)
  }

  if ("taxonomy_virus" %in% tables) {
    cat("  - Loading virus taxonomy...\n")
    virion_data$taxonomy_virus <- readr::read_csv(file.path(data_path, "taxonomy_virus.csv"), show_col_types = FALSE)
  }

  if ("provenance" %in% tables) {
    cat("  - Loading provenance data...\n")
    virion_data$provenance <- read_virion_csv(file.path(data_path, "provenance.csv.gz"))
  }

  if ("detection" %in% tables) {
    cat("  - Loading detection methods...\n")
    virion_data$detection <- vroom::vroom(file.path(data_path, "detection.csv.gz"), show_col_types = FALSE)
  }

  if ("temporal" %in% tables) {
    cat("  - Loading temporal data...\n")
    virion_data$temporal <- vroom::vroom(file.path(data_path, "temporal.csv.gz"), show_col_types = FALSE)
  }

  cat("Data loading complete!\n")
  virion_data
}

get_virion_versions <- function() {
  if (!require(virionData, quietly = TRUE)) {
    stop("virionData package not installed")
  }

  versions <- virionData::list_deposit_versions()
  summary_info <- virionData::deposit_summary()

  list(
    versions = versions,
    summary = summary_info
  )
}
