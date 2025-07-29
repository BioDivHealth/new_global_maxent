library(pacman)
p_load(here, rgbif, taxize, raster, dismo, 
      doParallel, rJava, XML, rgbif, Hmisc, readr, 
      stringr, purrr, dplyr, tidyr, magrittr, tidyverse)

library(dplyr)

load(file = "scripts/functions/wrld_simpl2.R")
source("scripts/New_functions/get_synonyms.R")
options(iucn_redlist_key="tiB4fspZ5oyjmPYd88F5NqpNFxitdb4mfqu4")

# Helper function from 0_SpList.R -----------------------------------------
collapse_vals <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x)])
  paste(x, collapse = sep)
}

# Load data ----------------------------------------------------------------
who_virion_hosts_short = read_csv(here("data_artur", "WHO", "virion", "who_pathogens_virion_hosts_summary.csv"))

host_species = who_virion_hosts_short %>%
  dplyr::select(Host, HostGenus, HostFamily, HostFlagID) %>%
  mutate(
    Host = str_to_sentence(Host),
    HostGenus = str_to_sentence(HostGenus),
    HostFamily = str_to_sentence(HostFamily)
  ) %>%
  distinct()

host_species_list = sort(unique(host_species$Host))
cat("Found", length(host_species_list), "unique host species to standardize\n")

hosts_ids = who_virion_hosts_short %>% 
  dplyr::select(Host, HostTaxID) %>%
  mutate(
    Host = str_to_sentence(Host)
  ) %>%
  distinct() %>% 
  group_by(Host) %>%
  summarise(HostTaxID = paste(unique(HostTaxID), collapse = ";"), .groups = "drop")


# Main species standardization loop ---------------------------------------
species_list = list()

cat("Starting species standardization process...\n")
for (i in progress:length(host_species_list)) {
#for (i in 1:5) {  # Remove this line once testing is complete
    sp = host_species_list[i]
    cat("Processing species", i, "of", length(host_species_list), ":", sp, "\n")
    
    species_list[[i]] = retrieve_syns_new(sp,   # [Character] The species name from which to collect taxonomic information
                                            n_times=3,  # [Numeric] Number of times the search is repeated until a data is found,default value = 1
                                            Gbif=TRUE)
    species_list[[i]]$type = "host"
    species_list[[i]]$host_species = sp    
}

progress = length(species_list) + 1


# Create one-row-per-species dataframe -----------------------------------
cat("Creating standardized taxonomic dataframe...\n")
tax_df <- map_dfr(species_list, function(rec) {
  
  ## 1. Summarise TaxDat -----------------------------------------------------
  td <- rec$TaxDat
  td_summary <- if (is.null(td) || nrow(td) == 0) {
    tibble()                       # no extra columns to add
  } else {
    td %>% summarise(across(everything(), collapse_vals), .groups = "drop")
  }
  
  ## 2. Scalar + collapsed vectors ------------------------------------------
  tibble(
    Submitted_name = rec$Submitted_name,
    correct_name = rec$correct_name,  # add the correct name to the record
    type = rec$type,  # add the type of species (host)
    taxon_level = rec$taxon_level,
    host_species = rec$host_species,  # add the original host species name
    Spp_syn        = collapse_vals(rec$Spp_syn),
    IUCN_spp       = collapse_vals(rec$IUCN_spp)
  ) %>%
    bind_cols(td_summary)          # add the TaxDat summary columns
})

# Add consolidated taxonomic information ----------------------------------
#tax_df = read_csv(here("data_artur", "WHO", "who_host_species_standardized.csv"))

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
tax_df$Genus = sapply(strsplit(as.character(tax_df$correct_name), " "), `[`, 1)


tax_df_joined = tax_df %>% 
  left_join(hosts_ids, by = c("Submitted_name" = "Host")) %>% 
  rename(Host = Submitted_name) %>%
  relocate(HostTaxID, .after = Host) %>% 
  relocate(Genus, Family, Order, Class, Phylum, .after = taxon_level)
# Create output directory if it doesn't exist

output_dir <- here("data_artur","WHO","virion")
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the standardized host species data --------------------------------
output_file <- here(output_dir, "who_host_species_standardized.csv")
write_csv(tax_df_joined, output_file)

cat("Standardization complete!\n")
cat("Processed", nrow(tax_df), "host species records\n")
cat("Results saved to:", output_file, "\n")

# Summary statistics ------------------------------------------------------
cat("\n=== STANDARDIZATION SUMMARY ===\n")
cat("Total species processed:", length(host_species_list), "\n")
cat("Records with correct names:", sum(!is.na(tax_df$correct_name)), "\n")
cat("Success rate:", round(sum(!is.na(tax_df$correct_name)) / nrow(tax_df) * 100, 1), "%\n")

# Taxonomic breakdown
cat("\nTaxonomic breakdown:\n")
if (sum(!is.na(tax_df$Phylum)) > 0) {
  cat("Phyla represented:\n")
  print(table(tax_df$Phylum, useNA = "ifany"))
}

if (sum(!is.na(tax_df$Class)) > 0) {
  cat("\nClasses represented:\n")
  print(table(tax_df$Class, useNA = "ifany"))
}

# Species with issues
problematic_species <- tax_df %>% 
  filter(is.na(correct_name)) %>%
  select(Submitted_name, host_species)

if (nrow(problematic_species) > 0) {
  cat("\nSpecies that could not be standardized:\n")
  print(problematic_species)
}

cat("\nProcess completed successfully!\n")


coalesce(tmp$IUCN_Class,  
         tmp$ITIS_Class,  
         tmp$GBIF_Class)

# Visualisations -----------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(here)
library(magrittr)
host = read_csv(here("data_artur","WHO","virion","who_host_species_standardized.csv"))
host %<>% dplyr::select(Host, correct_name, Genus, Family, Order, Class, Phylum)

# ------------------------------------------------------------------------------
# TAXONOMY VISUALIZATIONS
# ------------------------------------------------------------------------------

cat("\n=== CREATING TAXONOMY VISUALIZATIONS ===\n")

# Clean data for visualization
host_clean <- host %>%
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
    title = "Host Species Distribution by Class",
    subtitle = paste("Total:", nrow(host_clean), "host species"),
    fill = "Taxonomic Class"
  )

# 2. ORDER DISTRIBUTION BAR CHART (Top 15)
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
    title = "Top 20 Host Orders by Species Count",
    subtitle = "Colored by taxonomic class",
    x = "Taxonomic Order",
    y = "Number of Species",
    fill = "Class"
  )

print(p2_order_bars)



# Save all plots
output_dir <- here("figures", "taxonomy_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(here(output_dir, "host_class_distribution.png"), p1_class_pie, 
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(here(output_dir, "host_orders_top15.png"), p2_order_bars, 
       width = 12, height = 8, dpi = 300, bg = "white")
