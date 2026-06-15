# ------------------------------------------------------------------------------
# 2_2_CLOVER_Host_Clean.R
# ------------------------------------------------------------------------------
# Purpose: Standardize host species names from CLOVER database using the same
#          approach as used for VIRION hosts in 4_Host_Species_Clean.R
#
# Input  : CLOVER host associations from 2_1_CLOVER.R
# Output : Standardized taxonomic information for CLOVER host species
# ------------------------------------------------------------------------------|

library(pacman)
p_load(here, rgbif, taxize, raster, dismo, 
      doParallel, rJava, XML, rgbif, Hmisc, readr, 
      stringr, purrr, dplyr, tidyr, magrittr, tidyverse)

library(dplyr)

source(here("scripts", "associations", "working_inputs.R"))
source(here(
  "scripts",
  "associations",
  "network_building",
  "helpers",
  "legacy_who_compatibility_helpers.R"
))

load(file = "scripts/functions/wrld_simpl2.R")
source("scripts/New_functions/get_synonyms.R")
iucn_redlist_key <- Sys.getenv("IUCN_REDLIST_KEY", unset = Sys.getenv("IUCN_API_KEY", unset = ""))
if (nzchar(iucn_redlist_key)) {
  options(iucn_redlist_key = iucn_redlist_key)
}

# Load data ----------------------------------------------------------------
clover_hosts_file <- file.path(who_clover_dir, "who_bacteria_clover_hosts.csv")
clover_hosts_data <- read_csv(clover_hosts_file)

host_species <- clover_hosts_data %>%
  dplyr::select(Host, HostGenus, HostFamily, HostTaxID) %>%
  mutate(
    Host = str_to_sentence(Host),
    HostGenus = str_to_sentence(HostGenus),
    HostFamily = str_to_sentence(HostFamily)
  ) %>%
  distinct()

host_species_list <- sort(unique(host_species$Host))
cat("Found", length(host_species_list), "unique host species to standardize\n")

hosts_ids <- clover_hosts_data %>% 
  dplyr::select(Host, HostTaxID) %>%
  mutate(
    Host = str_to_sentence(Host)
  ) %>%
  distinct() %>% 
  group_by(Host) %>%
  summarise(HostTaxID = paste(unique(HostTaxID), collapse = ";"), .groups = "drop")

# Main species standardization loop ---------------------------------------
species_list <- list()

cat("Starting species standardization process...\n")

for (i in 1:length(host_species_list)) {
    sp <- host_species_list[i]
    cat("Processing species", i, "of", length(host_species_list), ":", sp, "\n")
    
    species_list[[i]] <- retrieve_syns_new(sp,   # [Character] The species name from which to collect taxonomic information
                                       n_times=10,  # [Numeric] Number of times the search is repeated until a data is found,default value = 1
                                       Gbif=TRUE)
    species_list[[i]]$type <- "host"
    species_list[[i]]$host_species <- sp    
}


# Create one-row-per-species dataframe -----------------------------------
cat("Creating standardized taxonomic dataframe...\n")
tax_df <- map_dfr(species_list, function(rec) {
  
  ## 1. Summarise TaxDat -----------------------------------------------------
  td <- rec$TaxDat
  td_summary <- if (is.null(td) || nrow(td) == 0) {
    tibble()                       # no extra columns to add
  } else {
    td %>% summarise(across(everything(), legacy_who_collapse_vals), .groups = "drop")
  }
  
  ## 2. Scalar + collapsed vectors ------------------------------------------
  tibble(
    Submitted_name = rec$Submitted_name,
    correct_name = rec$correct_name,  # add the correct name to the record
    type = rec$type,  # add the type of species (host)
    taxon_level = rec$taxon_level,
    host_species = rec$host_species,  # add the original host species name
    Spp_syn        = legacy_who_collapse_vals(rec$Spp_syn),
    IUCN_spp       = legacy_who_collapse_vals(rec$IUCN_spp)
  ) %>%
    bind_cols(td_summary)          # add the TaxDat summary columns
})

# Add consolidated taxonomic information ----------------------------------
tax_df$Phylum <- coalesce(tax_df$IUCN_Phylum, 
                          tax_df$ITIS_Phylum, 
                          tax_df$GBIF_Phylum)
tax_df$Class  <- coalesce(tax_df$IUCN_Class,  
                          tax_df$ITIS_Class,  
                          tax_df$GBIF_Class)
tax_df$Order  <- coalesce(tax_df$IUCN_Order,  
                          tax_df$ITIS_Order,  
                          tax_df$GBIF_Order)
tax_df$Family <- coalesce(tax_df$IUCN_Family, 
                          tax_df$ITIS_Family, 
                          tax_df$GBIF_Family)

# To get genus, we take first word from correct_name
tax_df$Genus <- sapply(strsplit(as.character(tax_df$correct_name), " "), `[`, 1)

tax_df_joined <- tax_df %>% 
  left_join(hosts_ids, by = c("Submitted_name" = "Host")) %>% 
  rename(Host = Submitted_name) %>%
  relocate(HostTaxID, .after = Host) %>% 
  relocate(Genus, Family, Order, Class, Phylum, .after = taxon_level)

for (i in 1:nrow(tax_df_joined)){
  if(nchar(tax_df_joined$Spp_syn[i]) > 0){
  tax_df_joined$Spp_syn[i] =  clean_synonyms2(tax_df_joined$Spp_syn[i])}
}

# Create output directory if it doesn't exist
output_dir <- who_clover_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the standardized host species data --------------------------------
output_file <- file.path(output_dir, "clover_host_species_standardized.csv")
write_csv(tax_df_joined, output_file)

cat("Standardization complete!\n")
cat("Processed", nrow(tax_df), "host species records\n")
cat("Results saved to:", output_file, "\n")
# ------------------------------------------------------------------------------
# TAXONOMY VISUALIZATIONS
# ------------------------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(here)
library(magrittr)

cat("\n=== CREATING TAXONOMY VISUALIZATIONS ===\n")

# Clean data for visualization
host_clean <- tax_df_joined %>%
  filter(!is.na(correct_name)) %>%
  mutate(
    Class = str_to_title(coalesce(Class, "Unknown")),
    Order = str_to_title(coalesce(Order, "Unknown")),
    Family = str_to_title(coalesce(Family, "Unknown")),
    Phylum = str_to_title(coalesce(Phylum, "Unknown"))
  )

# 1. CLASS DISTRIBUTION PIE CHART
p1_class_pie <- host_clean %>%
  count(Class, sort = TRUE) %>%
  mutate(
    percentage = round(n/sum(n) * 100, 1),
    label = paste0(Class, "\n(", n, " species, ", percentage, "%)")
  ) %>%
  ggplot(aes(x = "", y = n, fill = Class)) +
  geom_col(width = 1, color = "white", size = 0.5) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "CLOVER Host Species Distribution by Class",
    subtitle = paste("Total:", nrow(host_clean), "host species"),
    fill = "Taxonomic Class"
  )

# 2. ORDER DISTRIBUTION BAR CHART (Top 20)
p2_order_bars <- host_clean %>%
  count(Order, Class, sort = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(Order = fct_reorder(Order, n)) %>%
  ggplot(aes(x = Order, y = n, fill = Class)) +
  geom_col(alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  ) +
  labs(
    title = "Top 20 CLOVER Host Orders by Species Count",
    subtitle = "Colored by taxonomic class",
    x = "Taxonomic Order",
    y = "Number of Species",
    fill = "Class"
  )



# Display plots
cat("Displaying taxonomy visualizations...\n")
print(p1_class_pie)
print(p2_order_bars) 
