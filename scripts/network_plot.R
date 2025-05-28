library(igraph)
library(ggraph)
library(tidyverse)

setwd("~/DRTO_metaG/pathway_hog_connectivity")
filtered_results <- read.table(file = "high_corr_disease.txt", sep = "\t")
#filtered_results <- read.table(file = "high_corr_healthy.txt", sep = "\t")
#filtered_results <- read.table(file = "high_corr_AH.txt", sep = "\t")

# Create a graph from highly positively correlation pairs
# edges <- sig_results %>%
#   #filter(Correlation > 0.8) %>% #use this if you want to only plot positive correlations
#   filter(abs(Correlation) > 0.2) %>% #use this to plot both positive and negative correlations
#   dplyr::select(from = HOG_ID, to = Pathway, weight = Correlation)

edges <- filtered_results %>%
  #filter(Correlation > 0.8) %>% #use this if you want to only plot positive correlations
  filter(abs(Correlation) > 0.79) %>% #use this to plot both positive and negative correlations
  dplyr::select(from = HOG_ID, to = Pathway, weight = Correlation)

edges$color_weight <- edges$weight
edges$weight <- abs(edges$weight)

hog_nodes <- unique(edges$from)
pathway_nodes <- unique(edges$to)

# Create a node dataframe with type labels
nodes <- data.frame(
  name = unique(c(hog_nodes, pathway_nodes)),
  type = ifelse(unique(c(hog_nodes, pathway_nodes)) %in% hog_nodes, "HOG", "Pathway")
)

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

#I was getting a bunch of dyads so the following code removes them
#identify connected components
comps <- components(g)
# Get membership and component sizes
comp_sizes <- table(comps$membership)
# Get component IDs that are exactly 2 nodes in size
dyad_ids <- as.numeric(names(comp_sizes[comp_sizes == 2]))
# Find which nodes belong to those dyads
dyad_nodes <- names(comps$membership[comps$membership %in% dyad_ids])
# Remove those nodes (and their edges) from the graph
g_clean <- delete_vertices(g, dyad_nodes)
length(dyad_nodes) #26


## only positive correlations (i.e. when one increases in abundance so does the other)
network <- ggraph(g_clean, layout = "fr") +
  geom_edge_link(aes(edge_color = weight), edge_alpha = 1, edge_width = 1, show.legend = TRUE) +
  geom_node_point(aes(color = type), size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_gradient2(low = "blue", mid = "purple", high = "red", midpoint = 0.88) +
  scale_color_manual(values = c("HOG" = "orange", "Pathway" = "forestgreen")) +
  theme_void() +
  labs(title = "Clustered HOG–Pathway Correlation Network", color = "Node Type")
network 

ggsave("hog_pathway_network_dis_only_positive.svg", plot = network, width = 10, height = 8, units = "in", dpi = 300)

## all correlated clusters
network2 <- ggraph(g_clean, layout = "fr") +
  geom_edge_link(aes(edge_alpha = color_weight, edge_color = color_weight), show.legend = TRUE) +
  geom_node_point(aes(color = type), size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_color_manual(values = c("HOG" = "orange", "Pathway" = "forestgreen")) +
  theme_void() +
  labs(title = "Orthogroup–Pathway Correlation Network")
network2

ggsave("hog_pathway_network_dis.png", plot = network2, width = 10, height = 8, units = "in", dpi = 400)

## node metrics
node_metrics <- data.frame(
  Node = V(g)$name,
  Type = V(g)$type,
  Degree = degree(g),
  Betweenness = betweenness(g, normalized = TRUE),
  Closeness = closeness(g, normalized = TRUE)
)

# Optional: sort by degree descending
node_metrics <- node_metrics[order(-node_metrics$Degree), ]

# View first few rows
head(node_metrics)

# Average degree (mean number of edges per node)
avg_degree <- mean(degree(g))

# Graph density: how many edges exist vs how many possible
net_density <- edge_density(g)

# Modularity: requires community detection (e.g., fast greedy)
comm <- cluster_fast_greedy(g)
net_modularity <- modularity(comm)

# Combine into a summary
network_summary <- data.frame(
  Avg_Degree = avg_degree,
  Density = net_density,
  Modularity = net_modularity
)

network_summary

write.csv(node_metrics, "node_connectivity_metrics_dis.csv", row.names = FALSE)
write.csv(network_summary, "network_summary_metrics_dis.csv", row.names = FALSE)