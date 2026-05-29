library(pacman)
p_load(here, tidyverse, readr, magrittr)

source(here("scripts", "associations", "working_inputs.R"))

clover_network = read_csv(who_network_source_component_path("clover_who_network.csv"))
clover_disease_names = read_csv(who_network_source_component_path("clover_who_network.csv"))

virion_network = read_csv(who_network_source_component_path("virion_who_network.csv"))

clover_network$PathogenType = "bacteria"
virion_network$PathogenType = "virus"

clover_network <- clover_network %>%
  mutate(
    in_gibb_etal = if ("in_gibb_etal" %in% names(.)) in_gibb_etal else NA,
    in_empres_i = if ("in_empres_i" %in% names(.)) in_empres_i else NA,
    high_quality_detection = if ("high_quality_detection" %in% names(.)) high_quality_detection else DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing"),
    downstream_default_include = if ("downstream_default_include" %in% names(.)) downstream_default_include else high_quality_detection,
    downstream_review_reason = if ("downstream_review_reason" %in% names(.)) downstream_review_reason else NA_character_
  )

virion_network <- virion_network %>%
  mutate(
    high_quality_detection = if ("high_quality_detection" %in% names(.)) high_quality_detection else DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing"),
    downstream_default_include = if ("downstream_default_include" %in% names(.)) downstream_default_include else high_quality_detection,
    downstream_review_reason = if ("downstream_review_reason" %in% names(.)) downstream_review_reason else NA_character_
  )

names(clover_network)[which(!(names(clover_network) %in% names(virion_network)))]
names(virion_network)[which(!(names(virion_network) %in% names(clover_network)))]
clover_network %<>% select(-any_of("ID"))

combined_network = bind_rows(clover_network, virion_network)

combined_network %<>% rename(Host = Host_clean) %>%
  mutate(
    PathogenClass = tolower(PathogenClass),
    PathogenOrder = tolower(PathogenOrder),
    PathogenFamily = tolower(PathogenFamily),
    PathogenGenus = tolower(PathogenGenus),
    HostPhylum = tolower(HostPhylum),
    HostClass = tolower(HostClass),
    HostFamily = tolower(HostFamily),
    HostOrder = tolower(HostOrder)
    )

output_path <- who_raw_network_path()
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(combined_network, output_path)

unique_pairs = combined_network %>% select(Host, Pathogen) %>% distinct()
dim(unique_pairs)
names(combined_network)
# Preview random selection of hosts 

sample(unique(combined_network$Host), size  = 100)

domesticated = read_csv(who_network_domesticated_path())

which(unique(domesticated$scientific_name) %in% unique(combined_network$Host))
