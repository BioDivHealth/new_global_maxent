library(pacman)
p_load(here, tidyverse, readr, magrittr)

clover_network = read_csv(here("data_artur", "WHO", "networks", "clover_who_network.csv"))
virion_network = read_csv(here("data_artur", "WHO", "networks", "virion_who_network.csv"))

clover_network$PathogenType = "bacteria"
virion_network$PathogenType = "virus"

names(clover_network)[which(!(names(clover_network) %in% names(virion_network)))]
names(virion_network)[which(!(names(virion_network) %in% names(clover_network)))]
clover_network %<>% select(-ID)

combined_network = rbind(clover_network, virion_network)

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

unique_pairs = combined_network %>% select(Host, Pathogen) %>% distinct()
names(combined_network)