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

load(file = "scripts/functions/wrld_simpl2.R")
source("scripts/New_functions/get_synonyms.R")
options(iucn_redlist_key="tiB4fspZ5oyjmPYd88F5NqpNFxitdb4mfqu4")

# Helper function from 0_SpList.R -----------------------------------------
collapse_vals <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x)])
  paste(x, collapse = sep)
}

# Load data ----------------------------------------------------------------
clover_hosts_file <- here("data_artur", "WHO", "clover", "who_bacteria_clover_hosts.csv")
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
# Initialize progress variable (uncomment and set if resuming)
# progress <- 1

for (i in 1:length(host_species_list)) {
#for (i in 1:5) {  # Remove this line once testing is complete
    sp <- host_species_list[i]
    cat("Processing species", i, "of", length(host_species_list), ":", sp, "\n")
    
    species_list[[i]] <- retrieve_syns(sp,   # [Character] The species name from which to collect taxonomic information
                                       n_times=3,  # [Numeric] Number of times the search is repeated until a data is found,default value = 1
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

# Create output directory if it doesn't exist
output_dir <- here("data_artur", "WHO", "clover")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save the standardized host species data --------------------------------
output_file <- here(output_dir, "clover_host_species_standardized.csv")
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

# 3. HIERARCHICAL TREEMAP
if (!require(treemap)) install.packages("treemap")
library(treemap)

treemap_data <- host_clean %>%
  count(Class, Order, name = "species_count") %>%
  filter(species_count >= 2) %>%  # Only show orders with 2+ species
  arrange(desc(species_count))

if (nrow(treemap_data) > 0) {
  p3_treemap <- treemap(
    treemap_data,
    index = c("Class", "Order"),
    vSize = "species_count",
    type = "index",
    palette = "Set3",
    title = "CLOVER Host Taxonomic Diversity: Class > Order",
    fontsize.title = 14,
    fontsize.labels = c(12, 10),
    align.labels = list(c("left", "top"), c("center", "center")),
    overlap.labels = 0.5,
    inflate.labels = FALSE
  )
}

# 4. DETAILED TAXONOMIC BREAKDOWN
top_families <- host_clean %>%
  count(Class, Order, Family, sort = TRUE) %>%
  slice_head(n = 20)

if (nrow(top_families) > 0) {
  p4_family_detail <- top_families %>%
    mutate(
      Family = fct_reorder(Family, n),
      Class_Order = paste(Class, Order, sep = " - ")
    ) %>%
    ggplot(aes(x = Family, y = n, fill = Class)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    ) +
    labs(
      title = "Top 20 CLOVER Host Families by Species Count",
      subtitle = "Showing taxonomic class",
      x = "Family",
      y = "Number of Species",
      fill = "Class"
    )
}

# Save all plots
output_plots_dir <- here("figures", "clover_taxonomy_plots")
dir.create(output_plots_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(here(output_plots_dir, "clover_host_class_distribution.png"), p1_class_pie, 
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(here(output_plots_dir, "clover_host_orders_top20.png"), p2_order_bars, 
       width = 12, height = 8, dpi = 300, bg = "white")

if (exists("p4_family_detail")) {
  ggsave(here(output_plots_dir, "clover_host_families_top20.png"), p4_family_detail, 
         width = 12, height = 10, dpi = 300, bg = "white")
}

# Save treemap as PNG
if (nrow(treemap_data) > 0) {
  png(here(output_plots_dir, "clover_host_taxonomy_treemap.png"), width = 12, height = 8, 
      units = "in", res = 300, bg = "white")
  treemap(
    treemap_data,
    index = c("Class", "Order"),
    vSize = "species_count",
    type = "index",
    palette = "Set3",
    title = "CLOVER Host Taxonomic Diversity: Class > Order",
    fontsize.title = 14,
    fontsize.labels = c(12, 10),
    align.labels = list(c("left", "top"), c("center", "center")),
    overlap.labels = 0.5,
    inflate.labels = FALSE
  )
  dev.off()
}

# Display plots
cat("Displaying taxonomy visualizations...\n")
print(p1_class_pie)
print(p2_order_bars) 

if (exists("p4_family_detail")) {
  print(p4_family_detail)
}

# Taxonomic diversity summary
diversity_summary <- host_clean %>%
  group_by(Class) %>%
  summarise(
    Species = n(),
    Orders = n_distinct(Order),
    Families = n_distinct(Family),
    Genera = n_distinct(Genus),
    .groups = "drop"
  ) %>%
  arrange(desc(Species)) %>%
  mutate(
    `Avg Species/Order` = round(Species/Orders, 1),
    `Avg Species/Family` = round(Species/Families, 1)
  )

# Print summary table
cat("\n=== TAXONOMIC DIVERSITY SUMMARY ===\n")
print(diversity_summary)

if (nrow(top_families) > 0) {
  cat("\n=== TOP FAMILIES BY SPECIES COUNT ===\n")
  print(top_families)
}

cat("\nTaxonomy visualizations completed!\n")
cat("All plots saved to:", output_plots_dir, "\n")

# Species richness by taxonomic level
richness_stats <- tibble(
  Level = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  Count = c(
    n_distinct(host_clean$Phylum),
    n_distinct(host_clean$Class), 
    n_distinct(host_clean$Order),
    n_distinct(host_clean$Family),
    n_distinct(host_clean$Genus),
    nrow(host_clean)
  )
)

cat("\n=== TAXONOMIC RICHNESS ACROSS LEVELS ===\n")
print(richness_stats)

cat("\nVisualization analysis complete!\n") 