# ------------------------------------------------------------------------------|
# 6_Network_Examples.R  
# ------------------------------------------------------------------------------|
# Purpose: Focused examples of different network visualization approaches
#          for pathogen-host associations
# ------------------------------------------------------------------------------|

library(pacman)
p_load(here, tidyverse, igraph, ggraph, visNetwork, plotly)

# Quick data loading (modify paths as needed)
host_data <- read_csv(here("data_artur", "WHO","virion", "who_pathogens_virion_hosts_summary.csv"))

# Note: If you want to use standardized taxonomy, you'll need to handle case differences:
host_taxonomy <- read_csv(here("data_artur", "WHO","virion", "who_host_species_standardized.csv"))
host_data <- host_data %>%
  mutate(Host_lower = str_to_lower(Host)) %>%
  left_join(host_taxonomy %>%
    mutate(Host_lower = str_to_lower(Host)) %>%
    select(Host_lower, correct_name, Class),
    by = "Host_lower") %>%
  mutate(Host = coalesce(correct_name, Host))

# Prepare simple network data
simple_network <- host_data %>%
  select(Virus, Host, VirusFamily, HostClass, `PHEIC risk`) %>%
  filter(!is.na(Virus), !is.na(Host)) %>%
  # Simplify for example
  #slice_head(n = 100) %>%
  distinct()

# ------------------------------------------------------------------------------
# EXAMPLE 1: Basic Force-Directed Network
# ------------------------------------------------------------------------------
create_basic_network <- function(data) {
  # Create edge list
  edges <- data %>%
    count(Virus, Host, name = "weight") %>%
    select(from = Virus, to = Host, weight)
  
  # Create node list
  nodes <- tibble(
    name = unique(c(edges$from, edges$to)),
    type = ifelse(name %in% edges$from, "Pathogen", "Host")
  )
  
  # Create igraph object
  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  # Basic plot
  plot(g, 
       vertex.color = ifelse(V(g)$type == "Pathogen", "red", "lightblue"),
       vertex.shape = ifelse(V(g)$type == "Pathogen", "circle", "square"),
       vertex.size = 8,
       vertex.label.cex = 0.7,
       edge.width = E(g)$weight * 0.5,
       layout = layout_with_fr(g),
       main = "Basic Pathogen-Host Network")
  
  return(g)
}

# ------------------------------------------------------------------------------  
# EXAMPLE 2: Interactive visNetwork
# ------------------------------------------------------------------------------
create_interactive_example <- function(data) {
  # Prepare edges
  edges <- data %>%
    count(Virus, Host, name = "value") %>%
    rename(from = Virus, to = Host)
  
  # Prepare nodes
  all_names <- unique(c(edges$from, edges$to))
  nodes <- tibble(
    id = all_names,
    label = all_names,
    group = ifelse(id %in% edges$from, "Pathogen", "Host"),
    color = ifelse(group == "Pathogen", "#d62728", "#2ca02c"),
    size = 20
  )
  
  # Create interactive network
  visNetwork(nodes, edges) %>%
    visNodes(borderWidth = 2) %>%
    visEdges(arrows = "to") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = FALSE)
}

# ------------------------------------------------------------------------------
# EXAMPLE 3: Risk-Based Network (ggraph)
# ------------------------------------------------------------------------------
create_risk_network <- function(data) {
  # Filter to high-risk pathogens
  high_risk_data <- data %>%
    filter(str_detect(`PHEIC risk`, "High|Medium")) %>%
    slice_head(n = 50)  # Limit for clarity
  
  # Create network
  edges <- high_risk_data %>%
    select(from = Virus, to = Host, risk = `PHEIC risk`) %>%
    distinct()
  
  nodes <- tibble(
    name = unique(c(edges$from, edges$to)),
    type = ifelse(name %in% edges$from, "Pathogen", "Host")
  )
  
  # Create tidygraph object
  library(tidygraph)
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  
  # Plot with ggraph
  graph %>%
    ggraph(layout = "stress") +
    geom_edge_link(color = "grey70", alpha = 0.7) +
    geom_node_point(aes(color = type, size = type), alpha = 0.8) +
    scale_color_manual(values = c("Pathogen" = "#d62728", "Host" = "#2ca02c")) +
    scale_size_manual(values = c("Pathogen" = 4, "Host" = 3)) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE, max.overlaps = 10) +
    theme_void() +
    theme(legend.position = "bottom") +
    labs(title = "High-Risk Pathogen Networks",
         subtitle = "WHO priority pathogens and their hosts")
}

# ------------------------------------------------------------------------------
# EXAMPLE 4: Bipartite Network Layout
# ------------------------------------------------------------------------------
create_bipartite_network <- function(data) {
  # Select subset for clarity
  subset_data <- data %>%
    group_by(VirusFamily) %>%
    slice_head(n = 3) %>%
    ungroup() %>%
    slice_head(n = 30)
  
  # Create bipartite graph
  edges <- subset_data %>%
    select(Virus, Host) %>%
    distinct()
  
  # Create igraph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Set bipartite attribute
  V(g)$type <- V(g)$name %in% edges$Virus
  
  # Bipartite layout
  layout_matrix <- layout_as_bipartite(g)
  
  plot(g,
       layout = layout_matrix,
       vertex.color = ifelse(V(g)$type, "#ff7f0e", "#1f77b4"),
       vertex.shape = ifelse(V(g)$type, "circle", "square"),
       vertex.size = 8,
       vertex.label.cex = 0.7,
       main = "Bipartite Pathogen-Host Network")
}

# ------------------------------------------------------------------------------
# EXAMPLE 5: Centrality Analysis
# ------------------------------------------------------------------------------
analyze_centrality <- function(data) {
  # Create network
  edges <- data %>%
    count(Virus, Host) %>%
    select(from = Virus, to = Host, weight = n)
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Calculate centrality measures
  centrality_df <- tibble(
    name = V(g)$name,
    degree = degree(g),
    betweenness = betweenness(g),
    closeness = closeness(g),
    eigenvector = eigen_centrality(g)$vector
  ) %>%
    mutate(
      type = ifelse(name %in% edges$from, "Pathogen", "Host")
    ) %>%
    arrange(desc(betweenness))
  
  # Plot top central nodes
  top_nodes <- centrality_df %>%
    slice_head(n = 20)
  
  p <- ggplot(top_nodes, aes(x = reorder(name, betweenness), y = betweenness, fill = type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Pathogen" = "#d62728", "Host" = "#2ca02c")) +
    labs(title = "Most Central Nodes in Pathogen-Host Network",
         subtitle = "Ranked by betweenness centrality",
         x = "Node", y = "Betweenness Centrality") +
    theme_minimal()
  
  print(p)
  return(centrality_df)
}

# ------------------------------------------------------------------------------
# RUN EXAMPLES
# ------------------------------------------------------------------------------

cat("Creating network visualization examples...\n")

# Example 1: Basic network
cat("1. Creating basic network...\n")
basic_net <- create_basic_network(simple_network)

# Example 2: Interactive network (uncomment to run)
# cat("2. Creating interactive network...\n")
# interactive_net <- create_interactive_example(simple_network)
# interactive_net

# Example 3: Risk-based network
cat("3. Creating risk-based network...\n")
p_risk <- create_risk_network(simple_network)
print(p_risk)

# Example 4: Bipartite layout  
cat("4. Creating bipartite network...\n")
create_bipartite_network(simple_network)

# Example 5: Centrality analysis
cat("5. Analyzing network centrality...\n")
centrality_results <- analyze_centrality(simple_network)

cat("Examples completed! Try running individual functions with your full dataset.\n")

# ------------------------------------------------------------------------------
# TIPS FOR FURTHER EXPLORATION
# ------------------------------------------------------------------------------
cat("\n=== TIPS FOR NETWORK ANALYSIS ===\n")
cat("1. Filter by virus families of interest for focused networks\n")
cat("2. Use host taxonomic classes to group and color nodes\n") 
cat("3. Weight edges by detection confidence or number of studies\n")
cat("4. Create separate networks for different risk levels\n")
cat("5. Use community detection algorithms to find clusters\n")
cat("6. Export to Gephi or Cytoscape for advanced analysis\n") 