library(pacman)
p_load(here, tidyverse, readr, magrittr)

clover_network = read_csv(here("pathogen_association_data", "WHO",
                               "networks", "clover_who_network.csv"))
clover_disease_names = read_csv(here("pathogen_association_data", 
                                     "WHO", "networks", "clover_who_network.csv"))

virion_network = read_csv(here("pathogen_association_data",
                               "WHO", "networks", "virion_who_network.csv"))

clover_network$PathogenType = "bacteria"
virion_network$PathogenType = "virus"

clover_network <- clover_network %>%
  mutate(
    in_gibb_etal = if ("in_gibb_etal" %in% names(.)) in_gibb_etal else NA,
    in_empres_i = if ("in_empres_i" %in% names(.)) in_empres_i else NA
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

write_csv(combined_network, here("pathogen_association_data", "WHO", "networks", "combined_who_network.csv"))

unique_pairs = combined_network %>% select(Host, Pathogen) %>% distinct()
dim(unique_pairs)
names(combined_network)
# Preview random selection of hosts 

sample(unique(combined_network$Host), size  = 100)

domesticated = read_csv(here("pathogen_association_data","WHO","domesticated","domesticated_lab_farmed.csv"))

which(unique(domesticated$scientific_name) %in% unique(combined_network$Host))
