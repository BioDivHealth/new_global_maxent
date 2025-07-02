# Install required libraries ----
library(spocc)
library(here)
library(rgbif)
library(taxize)
library(raster)
library(dismo)
library(rgdal)
library(maptools)
library(doParallel)
library(rgeos)
library(rJava)
library(XML)
library(rgbif)
library(rgdal)
library(Hmisc)
library(readr)
library(stringr)
library(purrr)

# Read in data ----
# Make sure to read in with Latin1 encoding
spec1 <- read.csv(here("data", "disease_table28.csv"), stringsAsFactors = FALSE, 
                  fileEncoding = "Latin1")
options(iucn_redlist_key="tiB4fspZ5oyjmPYd88F5NqpNFxitdb4mfqu4")
spec1b <- spec1

# Load world map and helper functions
load(file = "scripts/functions/wrld_simpl2.R")
source("scripts/New_functions/get_synonyms.R")

## helper: paste unique, drop NA
collapse_vals <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x)])
  paste(x, collapse = sep)
}


for (i in 1:length(spec1b$name)){
  r = spec1b[i,]
  cat("Disease name:", r$name, "\n")
  cat("Standardised disease name:", r$disease, "\n")
  
  countries <- r$countries
  subs2 <- wrld_simpl2[wrld_simpl2$ISO2 %in% strsplit(countries, ",")[[1]], ]
  if (is.na(countries)) {
    next
  }
  
  dis <- r$name
  list_for_disease = list()
  #names(list_for_disease) = dis
  hosts <- capitalize(strsplit(r$all_host, ",")[[1]])
  vectors <- capitalize(strsplit(r$all_vectors, ",")[[1]])
  second <- strsplit(r$secondary, ",")[[1]]
  hv <- data.frame(species = c(hosts, vectors, second), 
                   type = c(rep("hosts", length(hosts)), 
                            rep("vectors", length(vectors)), 
                            rep("second", length(second))), 
                   stringsAsFactors = FALSE)
  hv <- hv[!is.na(hv$species), ]
  hv <- hv[!(hv$species == "NA"), ]
  dom <- hv[hv$species %in% c(c("human", "cattle", "ducks", "chickens", "pigs", "sheep", "goats"), 
                              capitalize(c("human", "cattle", "ducks", "chickens", "pigs", "sheep", "goats"))), ]
  hv2 <- hv[!hv$species %in% c(c("human", "cattle", "ducks", "chickens", "pigs", "sheep", "goats"), 
                               capitalize(c("human", "cattle", "ducks", "chickens", "pigs", "sheep", "goats"))), ]
  print(hv2)
  
  if (!nrow(dom) == 0) {
    domesticated = dom 
    domesticated$disease_og = r$name  # add the submitted name to the record
    domesticated$disease = r$disease  # add the submitted name to the record
    domesticated$Submitted_name = dom$species  # add the submitted name to the record
    domesticated$domesticated = TRUE
    domesticated = domesticated %>% select(-species)  # remove the species column, as it is already in Submitted_name
  } else {
    domesticated = data.frame()  # create an empty dataframe if no domesticated species found
  }
  
  if (!nrow(hv2) == 0) {
    cat("Species list for disease:", dis, ":\n", paste(hv2$species, collapse = ", "), "\n")
    for (j in 1:nrow(hv2)){
      
      sp = hv2$species[j]
      print(sp)
      list_for_disease[[j]] = retrieve_syns(sp,   # [Character] The species name from which to collect taxonomic information
                                            n_times=3,  # [Numeric] Number of times the search is repeated until a data is found,default value = 1
                                            Gbif=FALSE)
      list_for_disease[[j]]$type = hv2$type[j]  # add the type of species (host/vector/second) to the record
      list_for_disease[[j]]$disease_og = r$name  # add the submitted name to the record
      list_for_disease[[j]]$disease = r$disease  # add the submitted name to the record
    }
    ## one-row-per-species dataframe ──────────────────────────────────────────────
    tax_df <- map_dfr(list_for_disease, function(rec) {
      
      ## 1.  summarise TaxDat  -----------------------------------------------------
      td <- rec$TaxDat
      td_summary <- if (is.null(td) || nrow(td) == 0) {
        tibble()                       # no extra columns to add
      } else {
        td %>% summarise(across(everything(), collapse_vals), .groups = "drop")
      }
      
      ## 2.  scalar + collapsed vectors  ------------------------------------------
      tibble(
        Submitted_name = rec$Submitted_name,
        correct_name = rec$correct_name,  # add the correct name to the record
        type = rec$type,  # add the type of species (host/vector/second)
        taxon_level = rec$taxon_level,
        #ITIS_species_in_genus = rec$ITIS_species_in_genus,
        disease_og = rec$disease_og,  # add the submitted name to the record
        disease = rec$disease,  # add the submitted name to the record
        Spp_syn        = collapse_vals(rec$Spp_syn),
        IUCN_spp       = collapse_vals(rec$IUCN_spp)
      ) %>%
        bind_cols(td_summary)          # add the TaxDat summary columns
    })
  } else {
    tax_df <- data.frame()  # create an empty dataframe if no species found
  }
  
  print("#############################")
  print(i)
  print("#############################")
  
  combined_df <- bind_rows(tax_df, domesticated)
  
  # Save the combined dataframe to a CSV file
  disease_string = gsub(" ", "_", r$name)
  cat("Saving data for disease:", disease_string, "\n")
  output_file <- here("data", "diseases_species", paste0(disease_string, "_species.csv"))
  write.csv(combined_df, output_file, row.names = FALSE)
}


