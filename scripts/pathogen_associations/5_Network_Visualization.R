# ------------------------------------------------------------------------------|
# 5_Network_Visualization.R
# ------------------------------------------------------------------------------|
# Purpose: Create network visualizations of WHO pathogen-host associations
#          using standardized data from previous processing steps
#
# Input:   who_pathogens_virion_hosts_summary.csv (from 3_WHO_Virion_Hosts.R)
#          who_host_species_standardized.csv (from 4_Host_Species_Clean.R)
#
# Output:  Interactive and static network plots
# ------------------------------------------------------------------------------|

# ------------------------------| Load libraries |------------------------------
library(pacman)
p_load(here, tidyverse, igraph, ggraph, networkD3, visNetwork, 
       plotly, RColorBrewer, viridis, cowplot, scales)

# Additional network packages
if (!require(tidygraph)) install.packages("tidygraph")
library(tidygraph)

# ------------------------------| Load data |--------------------------------
cat("Loading pathogen-host association data...\n")

# Load the main association data
host_associations <- read_csv(here("data_artur", "WHO", "virion", "who_pathogens_virion_hosts_summary.csv"))
host_associations_long <- read_csv(here("data_artur", "WHO", "virion", "who_pathogens_virion_hosts_long.csv"))

# Harmonize Virus names to standardized taxonomy
synonyms <- c(
  "influenza a virus"                       = "alphainfluenzavirus influenzae",
  "chikungunya virus"                       = "alphavirus chikungunya",
  "dengue virus"                            = "orthoflavivirus denguei",
  "zika virus"                              = "orthoflavivirus zikaense",
  "west nile virus"                         = "orthoflavivirus nilense",
  "yellow fever virus"                      = "orthoflavivirus flavi",
  "monkeypox virus"                         = "orthopoxvirus monkeypox",
  "henipavirus nipahense"                   = "henipavirus nipahense",
  "middle east respiratory syndrome-related coronavirus" =
    "betacoronavirus cameli"
  # Leave SARS rows as they are unless you deliberately merge them
)
 host_associations$Virus_std <- unname(synonyms[ tolower(host_associations$Virus) ])
 host_associations$Virus_std[ is.na(host_associations$Virus_std) ] <- host_associations$Virus[ is.na(host_associations$Virus_std) ]

# Rename columns
host_associations$Virus_og = host_associations$Virus  # Keep original names for reference
host_associations$Virus = host_associations$Virus_std  # Use standardized names for analysis

# Load standardized host taxonomy  
host_taxonomy <- read_csv(here("data_artur", "WHO", "virion", "who_host_species_standardized.csv"))
host_taxonomy$Host_lower = str_to_lower(host_taxonomy$Host)

# Clean and prepare data for network analysis
network_data <- host_associations %>%
  # Create lowercase host names for matching
  mutate(Host_lower = str_to_lower(Host)) %>%
  # Use standardized host names if available (match on lowercase)
  left_join(host_taxonomy %>% select(Host, Host_lower, correct_name, Class, Family, Order), 
            by = "Host_lower") %>%
  mutate(
    # Use standardized name if available, otherwise original
    Host_clean = coalesce(correct_name, Host.x),  # Host.x is from host_associations
    # Simplify virus names for better visualization
    Virus_clean = str_remove(Virus, " sp\\.$|strain.*$"),
    # Create risk categories
    Risk_category = case_when(
      str_detect(`PHEIC risk`, "High") ~ "High Risk",
      str_detect(`PHEIC risk`, "Medium") ~ "Medium Risk", 
      str_detect(`PHEIC risk`, "Low") ~ "Low Risk",
      TRUE ~ "Unknown Risk"
    )
  ) %>%
  # Filter for high-quality detections
  #filter(DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing")) %>%
  # Remove uncertain host identifications if desired
  filter(!HostFlagID | is.na(HostFlagID)) %>%
  select(Virus_clean, Host_clean, VirusFamily, HostClass = Class, HostFamily = Family, HostOrder = Order,
         Risk_category) %>%
  distinct() %>%
  filter(!is.na(Virus_clean), !is.na(Host_clean))

cat("Prepared", nrow(network_data), "pathogen-host associations for visualization\n")


# ------------------------------| Create network objects |-------------------

#' Create network graph from pathogen-host data
#' @param data Data frame with pathogen-host associations
#' @param min_connections Minimum connections for inclusion
#' @return tidygraph object
create_pathogen_network <- function(data, min_connections = 1) {
  
  # Create edge list
  edges <- data %>%
    group_by(Virus_clean, Host_clean) %>%
    summarise(
      weight = n(),
      virus_family = first(VirusFamily),
      host_order = first(HostOrder),
      risk_category = first(Risk_category),
      .groups = "drop"
    ) %>%
    filter(weight >= min_connections) %>%
    rename(from = Virus_clean, to = Host_clean)
  
  # Create node attributes
  virus_nodes <- edges %>%
    select(name = from, virus_family, risk_category) %>%
    distinct() %>%
    mutate(
      type = "Pathogen",
      group = virus_family,
      risk = risk_category,
      size_metric = "virus_connections"
    )
  
  host_nodes <- edges %>%
    select(name = to, host_order) %>%
    distinct() %>%
    mutate(
      type = "Host",
      group = host_order,
      risk = "Host",
      size_metric = "host_connections"
    )
  
  nodes <- bind_rows(virus_nodes, host_nodes) %>%
    # Calculate degree centrality
    left_join(
      bind_rows(
        edges %>% count(from, name = "degree") %>% rename(name = from),
        edges %>% count(to, name = "degree") %>% rename(name = to)
      ) %>%
      group_by(name) %>%
      summarise(degree = sum(degree), .groups = "drop"),
      by = "name"
    ) %>%
    mutate(
      degree = ifelse(is.na(degree), 0, degree),
      # Size nodes by degree centrality
      size = scales::rescale(log1p(degree), to = c(3, 12))
    ) %>%
    # Collapse rarely occurring groups to 'Other' to avoid huge legends
    group_by(group) %>%
    mutate(group_plot = ifelse(n() < 20, "Other", group)) %>%
    ungroup()
  
  # Create tidygraph object
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  
  return(graph)
}

# ------------------------------| Static network visualizations |-------------

#' Create ggraph network plot
#' @param graph tidygraph object
#' @param layout Layout algorithm
#' @param title Plot title
plot_network_ggraph <- function(graph, layout = "stress", title = "Pathogen-Host Network") {
  
  # Color palettes
  virus_colors <- RColorBrewer::brewer.pal(min(8, graph %>% 
    activate(nodes) %>% 
    filter(type == "Pathogen") %>% 
    pull(group) %>% 
    n_distinct()), "Set2")
  
  host_colors <- RColorBrewer::brewer.pal(min(8, graph %>% 
    activate(nodes) %>% 
    filter(type == "Host") %>% 
    pull(group) %>% 
    n_distinct()), "Set1")
  
  # Identify top 5 hosts by degree (get both names and node IDs)
  top_hosts_info <- graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(type == "Host") %>%
    arrange(desc(degree)) %>%
    slice_head(n = 1) %>%
    mutate(node_id = row_number())
  
  top_hosts <- top_hosts_info$name
  top_host_ids <- which(graph %>% activate(nodes) %>% as_tibble() %>% pull(name) %in% top_hosts)

  p <- graph %>%
    # Add edge attributes for styling
    activate(edges) %>%
    mutate(
      is_top_host_edge = to %in% top_host_ids,  # Check 'to' against node IDs
      edge_color = ifelse(is_top_host_edge, "#2c7fb8", "grey70"),
      edge_alpha = ifelse(is_top_host_edge, 0.8, 0.3)
    ) %>%
    ggraph(layout = layout) +
    # Edges
    geom_edge_link(aes(alpha = edge_alpha, width = weight, color = edge_color), 
                   show.legend = FALSE) +
    scale_edge_color_identity() +
    scale_edge_width(range = c(0.1, 1.2)) +
    scale_edge_alpha_identity() +
    # Nodes
    geom_node_point(aes(size = size, color = group_plot, shape = type), 
                    alpha = 0.9) +
    scale_size_identity() +
    scale_shape_manual(values = c("Pathogen" = 16, "Host" = 17)) +
    scale_color_brewer(palette = "Set3", na.translate = FALSE) +
         # Labels for nodes
     # Pathogen labels (plain)
     geom_node_text(aes(label = ifelse(type == "Pathogen", name, "")), 
                    size = 3, fontface = "plain", repel = TRUE, max.overlaps = 20) +
     # Top host labels (bold, bigger)
     geom_node_text(aes(label = ifelse(name %in% top_hosts & type == "Host", name, "")), 
                    size = 4.5, fontface = "bold", color = "black", repel = TRUE, max.overlaps = 20) +
    # Theming
    theme_graph() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      legend.box = "horizontal"
    ) +
    guides(
      color = guide_legend(title = "Group", override.aes = list(size = 4)),
      shape = guide_legend(title = "Node Type", override.aes = list(size = 4))
    ) +
    labs(title = title)
  
  return(p)
}

#' Create risk-focused network visualization
plot_risk_network <- function(graph, title = "High-Risk Pathogen Networks") {
  
  # Filter to high-risk pathogens and their hosts
  high_risk_graph <- graph %>%
    activate(nodes) %>%
    filter(risk %in% c("High Risk", "Host")) %>%
    # Keep only connected components
    filter(degree > 0)
  
  # Risk color palette
  risk_colors <- c("High Risk" = "#d62728", "Host" = "#2ca02c")
  
  p <- high_risk_graph %>%
    ggraph(layout = "stress") +
    geom_edge_link(aes(alpha = weight), color = "grey50", show.legend = FALSE) +
    geom_node_point(aes(size = size, color = risk, shape = type), alpha = 0.8) +
    scale_color_manual(values = risk_colors) +
    scale_size_identity() +
    scale_shape_manual(values = c("Pathogen" = 16, "Host" = 17)) +
    geom_node_text(aes(label = ifelse(degree > 2, name, "")), 
                   size = 3, repel = TRUE, max.overlaps = 15) +
    theme_graph() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, color = "#d62728")
    ) +
    labs(title = title,
         subtitle = "Showing WHO high-risk pathogens and their known hosts")
  
  return(p)
}

# ------------------------------| Interactive visualizations |----------------

#' Create interactive network with visNetwork
#' @param graph tidygraph object
#' @param title Plot title
create_interactive_network <- function(graph, title = "Interactive Pathogen-Host Network") {
  
  # Prepare nodes for visNetwork
  nodes_vis <- graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(
      id = row_number(),
      label = name,
      title = paste0(
        "<b>", name, "</b><br>",
        "Type: ", type, "<br>",
        "Group: ", group, "<br>",
        "Connections: ", degree
      ),
      color = case_when(
        type == "Pathogen" & risk == "High Risk" ~ "#d62728",
        type == "Pathogen" & risk == "Medium Risk" ~ "#ff7f0e", 
        type == "Pathogen" & risk == "Low Risk" ~ "#2ca02c",
        type == "Pathogen" ~ "#1f77b4",
        type == "Host" ~ "#9467bd"
      ),
      shape = ifelse(type == "Pathogen", "dot", "triangle")
    ) %>%
    select(id, label, title, color, shape, size)
  
  # Prepare edges for visNetwork
  edges_vis <- graph %>%
    activate(edges) %>%
    as_tibble() %>%
    mutate(
      from_id = match(from, nodes_vis$label),
      to_id = match(to, nodes_vis$label),
      title = paste0(
        "Connection: ", from, " → ", to, "<br>",
        "Weight: ", weight, "<br>",
      ),
      width = scales::rescale(weight, to = c(1, 5))
    ) %>%
    select(from = from_id, to = to_id, title, width)
  
  # Create interactive network
  visNetwork(nodes_vis, edges_vis, main = title) %>%
    visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE),
               nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(enabled = TRUE, stabilization = FALSE) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visLegend(position = "right")
}

# ------------------------------| Analysis functions |------------------------

#' Analyze network properties
#' @param graph tidygraph object
analyze_network_properties <- function(graph) {
  
  cat("\n=== NETWORK ANALYSIS ===\n")
  
  # Basic network statistics
  n_nodes <- vcount(graph)
  n_edges <- ecount(graph)
  n_pathogens <- graph %>% activate(nodes) %>% filter(type == "Pathogen") %>% nrow()
  n_hosts <- graph %>% activate(nodes) %>% filter(type == "Host") %>% nrow()
  
  cat("Network size:\n")
  cat("  - Total nodes:", n_nodes, "\n")
  cat("  - Total edges:", n_edges, "\n")
  cat("  - Pathogens:", n_pathogens, "\n")
  cat("  - Hosts:", n_hosts, "\n")
  
  # Centrality measures
  graph_analyzed <- graph %>%
    activate(nodes) %>%
    mutate(
      betweenness = centrality_betweenness(),
      closeness = centrality_closeness(),
      eigenvector = centrality_eigen()
    )
  
  # Top central nodes
  cat("\nMost central pathogens (by betweenness):\n")
  top_pathogens <- graph_analyzed %>%
    activate(nodes) %>%
    filter(type == "Pathogen") %>%
    arrange(desc(betweenness)) %>%
    slice_head(n = 5) %>%
    as_tibble()
  
  print(top_pathogens %>% select(name, degree, betweenness))
  
  cat("\nMost connected hosts:\n")
  top_hosts <- graph_analyzed %>%
    activate(nodes) %>%
    filter(type == "Host") %>%
    arrange(desc(degree)) %>%
    slice_head(n = 5) %>%
    as_tibble()
  
  print(top_hosts %>% select(name, group, degree))
  
  return(graph_analyzed)
}

# ------------------------------| Main execution |----------------------------

cat("Creating pathogen-host network...\n")

# Create main network
main_network <- create_pathogen_network(network_data, min_connections = 1)

# Analyze network properties
analyzed_network <- analyze_network_properties(main_network)

# Create visualizations
cat("\nGenerating network visualizations...\n")

# 1. Overall network
p1 <- plot_network_ggraph(main_network, layout = "fr", 
                         title = "WHO Pathogen-Host Association Network")

p1

# 2. Risk-focused network  
p2 <- plot_risk_network(main_network)

# 3. Family-level network (simplified)
family_network <- network_data %>%
  group_by(VirusFamily, HostClass) %>%
  summarise(connections = n(), .groups = "drop") %>%
  filter(connections >= 2) %>%
  rename(from = VirusFamily, to = HostClass, weight = connections) %>%
  create_pathogen_network()

p3 <- plot_network_ggraph(family_network, layout = "stress",
                         title = "Virus Family - Host Class Networks")

# Save static plots
output_dir <- here("figures", "network_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(here(output_dir, "pathogen_host_network.png"), p1, 
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(here(output_dir, "high_risk_network.png"), p2,
       width = 12, height = 8, dpi = 300, bg = "white")
ggsave(here(output_dir, "family_class_network.png"), p3,
       width = 10, height = 8, dpi = 300, bg = "white")

# Create interactive visualization
cat("Creating interactive network visualization...\n")
interactive_net <- create_interactive_network(main_network)

# Save interactive plot
htmlwidgets::saveWidget(interactive_net, 
                       here(output_dir, "interactive_pathogen_network.html"),
                       selfcontained = TRUE)

cat("\nNetwork visualizations completed!\n")
cat("Static plots saved to:", output_dir, "\n")
cat("Interactive plot: interactive_pathogen_network.html\n")

# ------------------------------| Summary statistics |------------------------

# Create summary tables
network_summary <- analyzed_network %>%
  activate(nodes) %>%
  as_tibble() %>%
  group_by(type, group) %>%
  summarise(
    count = n(),
    avg_degree = mean(degree),
    max_degree = max(degree),
    .groups = "drop"
  ) %>%
  arrange(type, desc(count))

cat("\nNetwork composition by group:\n")
print(network_summary)

# Risk distribution
risk_summary <- analyzed_network %>%
  activate(nodes) %>%
  filter(type == "Pathogen") %>%
  as_tibble() %>%
  count(risk, sort = TRUE)

cat("\nPathogen risk distribution:\n")
print(risk_summary)

cat("\nAnalysis complete! Check the figures/network_plots/ directory for outputs.\n") 
