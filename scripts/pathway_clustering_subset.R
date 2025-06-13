library(tidyverse)

setwd("~/DRTO_metaG/pathway_hog_connectivity")

# Load functional abundance data (samples as rows, pathways as columns)
func <- read.csv("abund_transpose.csv", check.names = FALSE)

# Set sample class names as rownames and remove the now redundant 'Class' column
rownames(func) <- func$Class
func <- func %>% select(-Class)

# Subset to only disease margin samples (sample names end with 'Dis')
func <- func[grepl("(Dis)$", rownames(func)), ]

# Subset to only naive samples
#func <- func[grepl("(N)$", rownames(func)), ]

# Subset to only apparently healthy samples, including AH samples on D colonies (sample names end with 'N' or 'AH')
#func <- func[grepl("(Dis-AH|AH)$", rownames(func)), ]

# Remove the first column if it contains "Unassigned" pathway label (double-check this assumption)
func <- func %>% select(-1)

# Compute Spearman correlation matrix between pathways (columns)
func_cor_matrix <- cor(func, method = "spearman", use = "pairwise.complete.obs")

# Convert correlation matrix to a distance matrix (only clusters positively correlated terms)
func_dist <- as.dist(1 - func_cor_matrix)

# Perform hierarchical clustering of pathways
hc <- hclust(func_dist, method = "average")

# Cut the dendrogram to form clusters (r > 0.9 → h = 0.1)
clusters <- cutree(hc, h = 0.2)

# Create a dataframe mapping pathways to their cluster ID
func_clusters <- data.frame(Pathway = names(clusters), Cluster = clusters)

# Transpose the abundance matrix so that pathways are rows and samples are columns
func_t_df <- as.data.frame(t(func))
func_t_df$Pathway <- rownames(func_t_df)

# Join cluster assignments to transposed functional data
func_t_df <- left_join(func_t_df, func_clusters, by = "Pathway")

# Average the abundance of pathways within each cluster, for each sample
# (This assumes averaging makes biological sense which it does for this rel abund data—alternatively, you could sum)
func_cluster_avg <- func_t_df %>%
  select(-Pathway) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

# Create a label for each cluster listing all pathways it contains
func_clusters_members <- func_clusters %>%
  group_by(Cluster) %>%
  summarise(Pathway = paste(Pathway, collapse = ", "), .groups = "drop")

# Add cluster member labels to the averaged abundance matrix
func_cluster_avg_labeled <- as.data.frame(func_cluster_avg) %>%
  left_join(func_clusters_members, by = "Cluster")

rm(func_clusters_members)
rm(func_cluster_avg)

# Use cluster names as rownames and remove redundant 'Cluster' column
rownames(func_cluster_avg_labeled) <- paste0("Cluster_", func_cluster_avg_labeled$Cluster)
func_cluster_avg_labeled <- func_cluster_avg_labeled %>% select(-Cluster)

write.table(func_cluster_avg_labeled, file = "functional_pathway_clusters_all.txt", sep = "\t")
