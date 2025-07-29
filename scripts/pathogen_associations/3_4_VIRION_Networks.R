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
library(magrittr)

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
  left_join(host_taxonomy %>% select(Host, Host_lower, correct_name,Phylum, Class, Family, Order), 
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
  select(Pathogen = Virus_clean, Host_clean, HostTaxID, PathogenTaxID = VirusTaxID, 
         PathogenGenus = VirusGenus, PathogenFamily = VirusFamily, 
         PathogenOrder = VirusOrder, PathogenClass = VirusClass, HostPhylum = Phylum,
         HostClass = Class, HostFamily = Family, HostOrder = Order, DetectionMethod,
         `PHEIC risk`) %>%
  distinct() %>%
  mutate(MainSource = "VIRION") %>%
  filter(!is.na(Pathogen), !is.na(Host_clean))

cat("Prepared", nrow(network_data), "pathogen-host associations for visualization\n")
dir.create(here("data_artur", "WHO", "networks"), showWarnings = FALSE)
write_csv(network_data, here("data_artur", "WHO", "networks", "virion_who_network.csv"))

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
      size = scales::rescale(sqrt(degree), to = c(2, 8))  # Use sqrt for better scaling
    ) %>%
    # Collapse rarely occurring groups to 'Other' to avoid huge legends
    group_by(group) %>%
    mutate(group_plot = ifelse(n() < 20, "Other", group)) %>%
    ungroup()
  
  # Create tidygraph object
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  return(graph)
}

# ------------------------------| Static network visualizations |-------------

#' Create ggraph network plot
#' @param graph tidygraph object
#' @param layout Layout algorithm
#' @param title Plot title
plot_network_ggraph <- function(graph, layout = "stress", title = "Pathogen-Host Network", pathogen_color = "#1f77b4") {
  
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
    slice_head(n = 7) %>%
    mutate(node_id = row_number())
  
  top_hosts <- top_hosts_info$name
  top_host_ids <- which(graph %>% activate(nodes) %>% as_tibble() %>% pull(name) %in% top_hosts)

  p <- graph %>%
    # Add edge attributes for styling
    activate(edges) %>%
    mutate(
      is_top_host_edge = to %in% top_host_ids,
      edge_color = ifelse(is_top_host_edge, "#e31a1c", "grey70"),  # lighter grey
      edge_alpha = ifelse(is_top_host_edge, 0.3, 0.2)  # more transparent
    ) %>%
    ggraph(layout = layout) +
    # Edges
    geom_edge_link(aes(alpha = edge_alpha, width = weight, color = edge_color), 
                   show.legend = FALSE) +
    scale_edge_color_identity() +
    scale_edge_width(range = c(0.05, 0.5)) +
    scale_edge_alpha_identity() +
    # Nodes
    geom_node_point(aes(size = size, color = group_plot, shape = type), 
                    alpha = 0.9) +
    # Overlay pathogen nodes with a fixed color
    geom_node_point(
      data = function(x) { dplyr::filter(x, type == "Pathogen") },
      aes(size = size),
      shape = 17,
      color = pathogen_color,
      alpha = 0.9
    ) +
    # Highlight top host nodes with custom color
    geom_node_point(
      data = function(x) { dplyr::filter(x, name %in% top_hosts & type == "Host") },
      aes(size = 5),
      shape = 16,
      color = "#e31a1c",  # red color for top hosts
      alpha = 0.2
    ) +
    scale_size_identity() +
    scale_shape_manual(values = c("Pathogen" = 17, "Host" = 16)) +
    scale_color_manual(values = c(
      palette_by_class <- c(
        # Birds - cool / natural tones
        "ACCIPITRIFORMES"   = "#7fc97f",  # green (raptors)
        "ANSERIFORMES"      = "#beaed4",  # lavender (ducks, geese)
        "CHARADRIIFORMES"   = "#a6d854",  # yellow-green (shorebirds)
        "GALLIFORMES"       = "#66c2a5",  # teal (chickens, pheasants)
        "PASSERIFORMES"     = "#b3cde3",  # light blue
        "PELECANIFORMES"    = "skyblue1",  # grey-blue
        
        # Mammals - warm tones
        "ARTIODACTYLA"      = "#FF665A",  # light orange (even-toed ungulates)
        "CARNIVORA"         = "#FF9FDF",  # red-orange
        "CHIROPTERA"        = "#FFBFB5",  # warm peach
        "PRIMATES"          = "#Ff473B",  # magenta
        "RODENTIA"          = "#FFC60C",  # deep orange
        
        # Other
        "Other"             = "#666666"   # dark grey
      )
    ), na.translate = FALSE) +
         # Labels for nodes
     # Pathogen labels (plain)
     geom_node_text(aes(label = ifelse(type == "Pathogen" & degree < 17, name, "")), 
                    size = 2, color = "black", repel = TRUE, max.overlaps = 10, alpha = 0.8) +
     geom_node_text(aes(label = ifelse(type == "Pathogen" & degree >= 17, name, "")), 
                    size = 2, fontface = "bold", color = "#1f77b4", repel = TRUE, max.overlaps = 30, alpha = 0.9) +
     # Top host labels (bold, bigger)
     geom_node_text(aes(label = ifelse(name %in% top_hosts & type == "Host", name, "")), 
                    size =4, fontface = "bold", color = "#e31a1c", repel = TRUE, max.overlaps = 20, alpha = 0.8) +
    # Theming
    theme_graph() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      legend.box = "vertical"
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
                   size = 2, repel = TRUE, max.overlaps = 20) +
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
  n_pathogens <- graph %>% activate(nodes) %>% as_tibble() %>% filter(type == "Pathogen") %>% nrow()
  n_hosts <- graph %>% activate(nodes) %>% as_tibble() %>% filter(type == "Host") %>% nrow()
  
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
#network_data %<>% filter(Host_clean!="Homo sapiens")
#network_data %<>% filter(Host_clean!="Macroglossus sobrinus")

main_network <- create_pathogen_network(network_data, min_connections = 1)

# Analyze network properties
analyzed_network <- analyze_network_properties(main_network)

# Create visualizations
cat("\nGenerating network visualizations...\n")

layouts = c("fr", "stress","auto", "kk")

for (type_l in layouts){
# 1. Overall network
p1 <- plot_network_ggraph(main_network, layout = type_l, 
                         title = "WHO Priority + VIRION Pathogen-Host Association Network")

p1

output_dir <- here("figures", "network_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(here(output_dir, paste0("pathogen_host_network_", type_l,".png")), p1, 
       width = 18, height = 12, dpi = 600, bg = "white", scale = 0.9)
}

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

# ------------------------------| Advanced Network Analysis |------------------
cat("\n=== ADVANCED NETWORK ANALYSIS ===\n")
cat("Performing centrality, modularity, and sampling bias analysis...\n")

#' Comprehensive network analysis including modularity and bridge species identification
#' @param graph tidygraph object
#' @return list containing analysis results and plots
comprehensive_network_analysis <- function(graph) {
  
  # Convert to igraph for advanced analysis
  g_igraph <- graph %>% as.igraph()
  
  # Calculate comprehensive centrality measures
  # For bipartite directed networks, we need to be more careful with centrality measures
  graph_enhanced <- graph %>%
    activate(nodes) %>%
    mutate(
      # Basic centrality measures that work well with directed bipartite networks
      # For bipartite networks: total degree = indegree + outdegree
      indegree_centrality = centrality_degree(mode = "in"), 
      outdegree_centrality = centrality_degree(mode = "out"),
      degree_centrality = indegree_centrality + outdegree_centrality,
      
      # For betweenness, convert to undirected temporarily for meaningful results
      betweenness_centrality = {
        # Convert to undirected igraph and calculate betweenness directly
        g_undirected <- as.undirected(as.igraph(.))
        igraph::betweenness(g_undirected, normalized = TRUE)
      },
      
      # Closeness - use harmonic version which handles disconnected components better
      closeness_centrality = centrality_harmonic(normalized = TRUE),
      
      # PageRank works well for directed networks
      page_rank = centrality_pagerank(),
      
      # Network position measures 
      constraint = node_constraint(),  # Burt's structural constraint
      
      # For bipartite networks, clustering is always 0, so calculate local neighborhood density instead
      local_neighborhood_density = map_local_dbl(.f = function(neighborhood, ...) {
        # Calculate density of neighborhood
        n_nodes <- vcount(neighborhood)
        n_edges <- ecount(neighborhood)
        if (n_nodes <= 1) return(0)
        max_edges <- n_nodes * (n_nodes - 1) / 2  # undirected
        return(n_edges / max_edges)
      })
    )
  
  # Detect communities using multiple algorithms
  cat("Detecting network communities/modules...\n")
  
  # Convert to undirected for community detection
  g_undirected <- as.undirected(g_igraph, mode = "collapse")
  
  # Multiple community detection algorithms
  communities_louvain <- cluster_louvain(g_undirected)
  communities_walktrap <- cluster_walktrap(g_undirected)
  communities_infomap <- cluster_infomap(g_undirected)
  
  # Add community membership to graph
  graph_enhanced <- graph_enhanced %>%
    activate(nodes) %>%
    mutate(
      community_louvain = membership(communities_louvain),
      community_walktrap = membership(communities_walktrap), 
      community_infomap = membership(communities_infomap),
      # Calculate modularity for each algorithm
      modularity_louvain = modularity(communities_louvain),
      modularity_walktrap = modularity(communities_walktrap),
      modularity_infomap = modularity(communities_infomap)
    )
  
  # Print modularity scores
  cat("Community detection results:\n")
  cat("  - Louvain modularity:", round(modularity(communities_louvain), 3), 
      "with", length(communities_louvain), "communities\n")
  cat("  - Walktrap modularity:", round(modularity(communities_walktrap), 3), 
      "with", length(communities_walktrap), "communities\n")
  cat("  - Infomap modularity:", round(modularity(communities_infomap), 3), 
      "with", length(communities_infomap), "communities\n")
  
  return(list(
    graph = graph_enhanced,
    communities = list(
      louvain = communities_louvain,
      walktrap = communities_walktrap,
      infomap = communities_infomap
    )
  ))
}

#' Identify bridge species (high betweenness centrality hosts)
#' @param graph_data Enhanced graph with centrality measures
identify_bridge_species <- function(graph_data) {
  
  bridge_analysis <- graph_data %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(type == "Host") %>%
    arrange(desc(betweenness_centrality)) %>%
    mutate(
      # Use indegree for hosts (number of pathogens they host)
      host_degree = indegree_centrality,
      # Standardize centrality measures for comparison
      betweenness_z = scale(betweenness_centrality)[,1],
      degree_z = scale(host_degree)[,1],
      # Identify bridge species (high betweenness relative to degree)
      bridge_score = betweenness_z - degree_z,
      is_bridge = bridge_score > 1  # 1 standard deviation above expected
    )
  
  cat("\nBridge species analysis (hosts with high betweenness centrality):\n")
  cat("Top 10 bridge species:\n")
  print(bridge_analysis %>% 
    select(name, group, host_degree, betweenness_centrality, bridge_score) %>%
    head(10))
  
  return(bridge_analysis)
}

#' Identify potentially oversampled species
#' @param graph_data Enhanced graph with centrality measures
identify_oversampled_species <- function(graph_data) {
  
  # Analysis for both pathogens and hosts
  sampling_analysis <- graph_data %>%
    activate(nodes) %>%
    as_tibble() %>%
    group_by(type) %>%
    mutate(
      # Calculate expected degree based on group size
      group_size = n(),
      # Use appropriate degree measure for each type
      relevant_degree = ifelse(type == "Host", indegree_centrality, outdegree_centrality),
      expected_degree = median(relevant_degree),
      # Sampling bias indicators
      degree_ratio = relevant_degree / expected_degree,
      degree_z_score = scale(relevant_degree)[,1],
      # Flag potentially oversampled (unusually high connectivity)
      potentially_oversampled = degree_z_score > 2  # 2 standard deviations above mean
    ) %>%
    ungroup()
  
  cat("\nPotential sampling bias analysis:\n")
  
  # Oversampled pathogens
  oversampled_pathogens <- sampling_analysis %>%
    filter(type == "Pathogen", potentially_oversampled) %>%
    arrange(desc(degree_centrality))
  
  cat("Potentially oversampled pathogens (unusually high host connectivity):\n")
  print(oversampled_pathogens %>% 
    select(name, group, relevant_degree, degree_ratio, degree_z_score) %>%
    head(10))
  
  # Oversampled hosts
  oversampled_hosts <- sampling_analysis %>%
    filter(type == "Host", potentially_oversampled) %>%
    arrange(desc(relevant_degree))
  
  cat("\nPotentially oversampled hosts (unusually high pathogen connectivity):\n")
  print(oversampled_hosts %>% 
    select(name, group, relevant_degree, degree_ratio, degree_z_score) %>%
    head(10))
  
  return(sampling_analysis)
}

#' Create modularity visualization
#' @param graph_data Enhanced graph with community assignments
#' @param communities Community detection results
plot_modularity_network <- function(graph_data, communities, algorithm = "louvain") {
  
  community_col <- paste0("community_", algorithm)
  
  # Create color palette for communities
  n_communities <- graph_data %>% 
    activate(nodes) %>% 
    pull(!!sym(community_col)) %>% 
    max()
  
  community_colors <- rainbow(n_communities)
  
  p <- graph_data %>%
    ggraph(layout = "stress") +
    geom_edge_link(alpha = 0.3, color = "grey70") +
    geom_node_point(aes(size = degree_centrality, 
                       color = factor(!!sym(community_col)),
                       shape = type), 
                   alpha = 0.8) +
    scale_color_manual(values = community_colors) +
    scale_size_continuous(range = c(1, 8), name = "Degree") +
    scale_shape_manual(values = c("Pathogen" = 17, "Host" = 16)) +
    geom_node_text(aes(label = ifelse(betweenness_centrality > quantile(betweenness_centrality, 0.9), 
                                     name, "")), 
                   size = 2, repel = TRUE, max.overlaps = 15) +
    theme_graph() +
    theme(legend.position = "right") +
    labs(title = paste("Network Communities -", str_to_title(algorithm), "Algorithm"),
         subtitle = paste("Modularity =", round(communities[[algorithm]]$modularity, 3)),
         color = "Community")
  
  return(p)
}

#' Create centrality comparison plot
#' @param graph_data Enhanced graph with centrality measures
plot_centrality_analysis <- function(graph_data) {
  
  centrality_data <- graph_data %>%
    activate(nodes) %>%
    as_tibble() %>%
    select(name, type, group, degree_centrality, betweenness_centrality, 
           closeness_centrality, page_rank) %>%
    pivot_longer(cols = c(ends_with("_centrality"), page_rank), 
                names_to = "centrality_type", 
                values_to = "centrality_value") %>%
    mutate(
      centrality_type = str_remove(centrality_type, "_centrality"),
      centrality_type = str_to_title(centrality_type)
    )
  
  # Centrality correlation plot
  p1 <- graph_data %>%
    activate(nodes) %>%
    as_tibble() %>%
    ggplot(aes(x = degree_centrality, y = betweenness_centrality)) +
    geom_point(aes(color = type, size = page_rank), alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("Pathogen" = "#1f77b4", "Host" = "#ff7f0e")) +
    labs(title = "Centrality Relationships",
         subtitle = "Degree vs Betweenness Centrality (Bridge Species Analysis)",
         x = "Degree Centrality", 
         y = "Betweenness Centrality",
         size = "PageRank\nScore") +
    theme_minimal()
  
  # Centrality distribution by type
  p2 <- centrality_data %>%
    ggplot(aes(x = centrality_type, y = centrality_value, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("Pathogen" = "#1f77b4", "Host" = "#ff7f0e")) +
    facet_wrap(~type, scales = "free_y") +
    labs(title = "Centrality Measure Distributions",
         x = "Centrality Measure", 
         y = "Centrality Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(correlation = p1, distribution = p2))
}

# XX Execute comprehensive analysis -------------
cat("Running comprehensive network analysis...\n")
network_data %<>% filter(Host_clean!="Macroglossus sobrinus")
main_network <- create_pathogen_network(network_data, min_connections = 1)
enhanced_results <- comprehensive_network_analysis(main_network)
enhanced_graph <- enhanced_results$graph
communities <- enhanced_results$communities

# Identify bridge species and oversampled species
bridge_species <- identify_bridge_species(enhanced_graph)
sampling_analysis <- identify_oversampled_species(enhanced_graph)

# Create advanced visualizations
cat("Creating advanced network visualizations...\n")

# Modularity plots for different algorithms
p_mod_louvain <- plot_modularity_network(enhanced_graph, communities, "louvain")
p_mod_walktrap <- plot_modularity_network(enhanced_graph, communities, "walktrap")  
p_mod_infomap <- plot_modularity_network(enhanced_graph, communities, "infomap")

# Centrality analysis plots
centrality_plots <- plot_centrality_analysis(enhanced_graph)

# Save advanced analysis plots
analysis_dir <- here("figures", "network_plots", "advanced_analysis")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(here(analysis_dir, "communities_louvain.png"), p_mod_louvain,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(here(analysis_dir, "communities_walktrap.png"), p_mod_walktrap,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(here(analysis_dir, "communities_infomap.png"), p_mod_infomap,
       width = 14, height = 10, dpi = 300, bg = "white")

ggsave(here(analysis_dir, "centrality_correlation.png"), centrality_plots$correlation,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(here(analysis_dir, "centrality_distributions.png"), centrality_plots$distribution,
       width = 12, height = 8, dpi = 300, bg = "white")

# Save analysis results to CSV files
results_dir <- here("data_artur", "WHO", "virion", "network_analysis")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Enhanced node metrics
enhanced_node_data <- enhanced_graph %>%
  activate(nodes) %>%
  as_tibble()

write_csv(enhanced_node_data, here(results_dir, "node_centrality_metrics.csv"))

# Bridge species results
write_csv(bridge_species, here(results_dir, "bridge_species_analysis.csv"))

# Sampling bias analysis
write_csv(sampling_analysis, here(results_dir, "sampling_bias_analysis.csv"))

# Community membership
community_membership <- enhanced_node_data %>%
  select(name, type, group, starts_with("community_"), starts_with("modularity_"))

write_csv(community_membership, here(results_dir, "community_membership.csv"))

# Summary statistics
cat("\n=== ADVANCED ANALYSIS SUMMARY ===\n")

# Key bridge species
top_bridges <- bridge_species %>%
  filter(is_bridge) %>%
  arrange(desc(bridge_score)) %>%
  head(5)

cat("Top 5 bridge species (critical for pathogen ecology connectivity):\n")
for(i in 1:nrow(top_bridges)) {
  cat(sprintf("  %d. %s (%s) - Bridge score: %.2f\n", 
              i, top_bridges$name[i], top_bridges$group[i], top_bridges$bridge_score[i]))
}

# Oversampling summary
oversampled_summary <- sampling_analysis %>%
  group_by(type) %>%
  summarise(
    total_species = n(),
    potentially_oversampled = sum(potentially_oversampled),
    proportion_oversampled = mean(potentially_oversampled),
    .groups = "drop"
  )

cat("\nSampling bias summary:\n")
print(oversampled_summary)

# Community summary
community_summary <- enhanced_node_data %>%
  group_by(type, community_louvain) %>%
  summarise(
    count = n(),
    avg_degree = mean(degree_centrality),
    .groups = "drop"
  ) %>%
  arrange(community_louvain, desc(count))

cat("\nCommunity composition (Louvain algorithm):\n")
print(community_summary)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Advanced network analysis files saved to:", results_dir, "\n")
cat("Advanced visualization plots saved to:", analysis_dir, "\n")
cat("\nKey findings:\n")
cat("- Bridge species help identify critical hosts for cross-pathogen transmission\n")
cat("- Community structure reveals pathogen-host interaction modules\n") 
cat("- Sampling bias analysis highlights potentially over-studied species\n")
cat("- Use these results to guide future sampling strategies and identify research gaps\n") 
