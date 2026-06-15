# -----------------------------------------------------------------------------|
# master_plus_host_network_helpers.R ----
# -----------------------------------------------------------------------------|
# Purpose: Shared low-level helpers for master-plus host-network stage scripts.
# -----------------------------------------------------------------------------|

host_network_clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

host_network_clean_key <- function(x) {
  key <- host_network_clean_text(x)
  key <- stringr::str_to_lower(key)
  key <- stringr::str_replace_all(key, "&", " and ")
  key <- stringr::str_replace_all(key, "[^a-z0-9]+", " ")
  stringr::str_squish(key)
}

host_network_split_semicolon_values <- function(x) {
  x <- host_network_clean_text(x)
  if (is.na(x)) {
    return(character(0))
  }

  values <- stringr::str_split(x, ";", simplify = FALSE)[[1]]
  values <- stringr::str_squish(values)
  values <- purrr::discard(values, ~ .x == "")
  unique(values)
}

host_network_first_non_missing <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    NA_character_
  } else {
    x[[1]]
  }
}

host_network_collapse_reasons <- function(...) {
  reasons <- c(...)
  reasons <- reasons[!is.na(reasons) & reasons != ""]
  if (length(reasons) == 0) {
    NA_character_
  } else {
    paste(unique(reasons), collapse = "; ")
  }
}

host_network_collapse_unique <- function(x) {
  x <- host_network_clean_text(x)
  x <- sort(unique(stats::na.omit(x)))

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

host_network_add_missing_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  for (col in missing) {
    data[[col]] <- NA
  }
  data
}

host_network_association_key <- function(data) {
  paste(
    host_network_clean_key(data$Disease_name),
    host_network_clean_text(data$PathogenTaxID),
    host_network_clean_key(data$Pathogen),
    host_network_clean_text(data$HostTaxID),
    host_network_clean_key(data$Host),
    sep = "|||"
  )
}

host_network_is_true <- function(x) {
  x %in% c(TRUE, "TRUE", "true", "True", 1, "1")
}

host_network_collapse_true_flag <- function(x) {
  values <- x[!is.na(x)]

  if (length(values) == 0) {
    return(NA)
  }

  any(host_network_is_true(values))
}

host_network_match_one_query <- function(row_df, source_links) {
  source_name <- row_df$host_query_source[[1]]
  source_subset <- source_links %>% dplyr::filter(source == source_name)
  if (nrow(source_subset) == 0) {
    return(tibble::tibble())
  }

  taxids <- row_df$query_taxids[[1]]
  pathogen_keys <- row_df$query_pathogen_keys[[1]]

  matched <- source_subset %>% dplyr::mutate(match_method = NA_character_)

  if (length(taxids) > 0) {
    matched <- matched %>%
      dplyr::filter(!is.na(source_pathogen_taxid), source_pathogen_taxid %in% taxids) %>%
      dplyr::mutate(match_method = "taxid")
  } else if (length(pathogen_keys) > 0) {
    matched <- matched %>%
      dplyr::filter(!is.na(source_pathogen_key), source_pathogen_key %in% pathogen_keys) %>%
      dplyr::mutate(match_method = "name")
  } else {
    matched <- tibble::tibble()
  }

  if (nrow(matched) == 0) {
    return(tibble::tibble())
  }

  matched %>%
    dplyr::mutate(
      analysis_unit_id = row_df$analysis_unit_id[[1]],
      master_row = row_df$master_row[[1]],
      disease_master_name = row_df$disease_master_name[[1]],
      resolved_disease_name = row_df$resolved_disease_name[[1]],
      resolved_pathogen_name = row_df$resolved_pathogen_name[[1]],
      resolved_pathogen_rank = row_df$resolved_pathogen_rank[[1]],
      preferred_match_source = row_df$preferred_match_source[[1]],
      host_query_bucket = row_df$host_query_bucket[[1]],
      host_query_include_default = row_df$host_query_include_default[[1]],
      host_query_source = row_df$host_query_source[[1]],
      host_query_pathogen_names = row_df$host_query_pathogen_names[[1]],
      host_query_taxids = row_df$host_query_taxids[[1]],
      match_review_flag = row_df$match_review_flag[[1]],
      shared_species_proxy_flag = row_df$shared_species_proxy_flag[[1]],
      match_review_notes = row_df$match_review_notes[[1]]
    )
}

host_network_read_host_standardization <- function(path, source_name) {
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    dplyr::mutate(
      dplyr::across(where(is.character), host_network_clean_text),
      host_source = source_name,
      HostTaxID = host_network_clean_text(HostTaxID),
      raw_host_key = host_network_clean_key(Host),
      std_host = dplyr::coalesce(host_network_clean_text(correct_name), host_network_clean_text(Host)),
      std_host_phylum = host_network_clean_text(Phylum),
      std_host_class = host_network_clean_text(Class),
      std_host_family = host_network_clean_text(Family),
      std_host_order = host_network_clean_text(Order)
    ) %>%
    dplyr::select(
      host_source,
      HostTaxID,
      raw_host_key,
      std_host,
      std_host_phylum,
      std_host_class,
      std_host_family,
      std_host_order
    ) %>%
    dplyr::filter(!is.na(std_host))
}

host_network_summarise_host_lookup <- function(data, group_cols, method_name, suffix) {
  clean_col <- paste0("clean_host_", suffix)
  phylum_col <- paste0("HostPhylum_", suffix)
  class_col <- paste0("HostClass_", suffix)
  family_col <- paste0("HostFamily_", suffix)
  order_col <- paste0("HostOrder_", suffix)

  data %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(group_cols), ~ !is.na(.x))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_std_hosts = dplyr::n_distinct(std_host),
      "{clean_col}" := host_network_first_non_missing(std_host),
      "{phylum_col}" := host_network_first_non_missing(std_host_phylum),
      "{class_col}" := host_network_first_non_missing(std_host_class),
      "{family_col}" := host_network_first_non_missing(std_host_family),
      "{order_col}" := host_network_first_non_missing(std_host_order),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_std_hosts == 1) %>%
    dplyr::select(-n_std_hosts) %>%
    dplyr::mutate("{paste0('method_', suffix)}" := method_name)
}

host_network_read_source_network <- function(path, host_network_source, source_table, required_common) {
  data <- readr::read_csv(path, show_col_types = FALSE, na = c("", "NA")) %>%
    dplyr::mutate(
      dplyr::across(where(is.character), host_network_clean_text),
      PathogenTaxID = host_network_clean_text(PathogenTaxID),
      HostTaxID = host_network_clean_text(HostTaxID),
      high_quality_detection = dplyr::coalesce(high_quality_detection, FALSE),
      downstream_default_include = dplyr::coalesce(downstream_default_include, FALSE),
      host_network_source = host_network_source,
      source_table = source_table
    )

  missing_required <- setdiff(required_common, names(data))
  if (length(missing_required) > 0) {
    stop(
      source_table,
      " missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  data <- host_network_add_missing_columns(
    data,
    c(
      "Host_raw",
      "host_name_cleaning_method",
      "source_database",
      "source_assoc_id",
      "source_host_flag_id",
      "host_taxonomy_ready",
      "host_taxonomy_flag",
      "is_human_host",
      "is_model_or_lab_host",
      "is_domestic_or_livestock_hint"
    )
  )

  data %>%
    dplyr::mutate(
      Host_raw = dplyr::coalesce(Host_raw, Host),
      host_name_cleaning_method = dplyr::coalesce(host_name_cleaning_method, "existing_who_network_host"),
      source_database = dplyr::coalesce(source_database, MainSource),
      host_taxonomy_flag = dplyr::case_when(
        !is.na(host_taxonomy_flag) ~ host_taxonomy_flag,
        is.na(Host) ~ "missing_name",
        is.na(HostTaxID) ~ "missing_taxid",
        stringr::str_detect(stringr::str_to_lower(Host), "\\b(sp|spp|species|unidentified|unknown|uncultured)\\b\\.?") ~ "unresolved_sp",
        stringr::str_detect(stringr::str_to_lower(Host), "^[a-z][a-z-]+\\s+[a-z][a-z.-]+(\\s+[a-z][a-z.-]+)?$") ~ "species_like",
        TRUE ~ "unresolved_sp"
      ),
      host_taxonomy_ready = dplyr::coalesce(
        host_taxonomy_ready,
        host_taxonomy_flag == "species_like" & !is.na(HostTaxID)
      ),
      is_human_host = dplyr::coalesce(
        is_human_host,
        stringr::str_to_lower(Host) == "homo sapiens" | HostTaxID == "9606"
      ),
      is_model_or_lab_host = dplyr::coalesce(
        is_model_or_lab_host,
        stringr::str_detect(
          stringr::str_to_lower(Host),
          paste(
            c(
              "^homo sapiens$",
              "^mus musculus$",
              "^rattus norvegicus$",
              "^rattus rattus$",
              "^cavia porcellus$",
              "^mesocricetus auratus$",
              "^oryctolagus cuniculus$",
              "^macaca\\b",
              "^chlorocebus\\b",
              "^callithrix\\b",
              "^gallus gallus$"
            ),
            collapse = "|"
          )
        )
      ),
      is_domestic_or_livestock_hint = dplyr::coalesce(
        is_domestic_or_livestock_hint,
        stringr::str_detect(
          stringr::str_to_lower(Host),
          paste(
            c(
              "^bos taurus$",
              "^bos indicus$",
              "^bubalus bubalis$",
              "^ovis aries$",
              "^capra hircus$",
              "^sus scrofa$",
              "^equus caballus$",
              "^equus asinus$",
              "^camelus\\b",
              "^lama glama$",
              "^alpaca$",
              "^vicugna pacos$",
              "^gallus gallus$",
              "^meleagris gallopavo$",
              "^anas platyrhynchos$",
              "^anas platyrhynchos domesticus$",
              "^canis lupus familiaris$",
              "^felis catus$"
            ),
            collapse = "|"
          )
        )
      )
    )
}

host_network_make_scope_aliases <- function(data, source_priority) {
  alias_cols <- intersect(
    c("analysis_unit", "analysis_unit_label", "source_pathogen", "source_previous_name", "source_msl39_viral_name"),
    names(data)
  )

  data %>%
    dplyr::select(
      dplyr::any_of(c("source_disease_name", "modelling_scope_status", "modelling_scope_reason")),
      dplyr::all_of(alias_cols)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(alias_cols),
      names_to = "pathogen_alias_source",
      values_to = "pathogen_alias"
    ) %>%
    dplyr::transmute(
      scope_priority = source_priority,
      disease_key = host_network_clean_key(source_disease_name),
      pathogen_key = host_network_clean_key(pathogen_alias),
      modelling_scope_status,
      modelling_scope_reason
    )
}
