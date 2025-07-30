library(pacman)
p_load(here, tidyverse, readr, magrittr)

clover_network = read_csv(here("pathogen_association_data", "WHO", "networks", "clover_who_network.csv"))
clover_disease_names = read_csv(here("pathogen_association_data", "WHO", "networks", "clover_who_network.csv"))

virion_network = read_csv(here("pathogen_association_data", "WHO", "networks", "virion_who_network.csv"))

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

combined_network$Pathogen = str_to_sentence(combined_network$Pathogen)
write_csv(combined_network, here("pathogen_association_data", "WHO", "networks", "combined_who_network.csv"))

unique_pairs = combined_network %>% select(Host, Pathogen) %>% distinct()
dim(unique_pairs)
names(combined_network)
# Preview random selection of hosts 

sample(unique(combined_network$Host), size  = 100)

domesticated = read_csv(here("pathogen_association_data","WHO","domesticated","domesticated_lab_farmed.csv"))

which(unique(domesticated$scientific_name) %in% unique(combined_network$Host))
