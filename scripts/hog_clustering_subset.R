library(tidyverse)

setwd("~/DRTO_metaG/pathway_hog_connectivity")
# Load Lauren's data (samples as rows, pathways as columns)
# ~7500 single copy orthologs conserved across all four species 
# Reads were imported into R using tximport, then normalized with DESeq2’s vst normalization method.
# I added a column to match her metadata with my class groups
HOGs = read.csv(file = "normalized_reads_corr_GK.csv", header = TRUE, )

# Average data within each class
averaged_HOGs <- HOGs %>%
  group_by(GK_class) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
averaged_HOGs <- as.data.frame(averaged_HOGs)
rm(HOGs)

# Set sample class names as rownames and remove the now redundant 'Class' column
rownames(averaged_HOGs) <- averaged_HOGs$GK_class
averaged_HOGs <- averaged_HOGs %>% select(-GK_class)

# Subset to only disease margin samples (sample names end with 'Dis')
#averaged_HOGs <- averaged_HOGs[grepl("Dis$", rownames(averaged_HOGs)), ]

#Subset to only naive samples (sample names end with 'N')
#averaged_HOGs <- averaged_HOGs[grepl("N$", rownames(averaged_HOGs)), ]

# Subset to only apparently healthy samples, including AH samples on D colonies (sample names end with 'N' or 'AH')
averaged_HOGs <- averaged_HOGs[grepl("(Dis-AH|AH)$", rownames(averaged_HOGs)), ]

# Compute Spearman correlation matrix between Hierarchical Orthologous Groups (HOGs) (columns)
go_cor_matrix <- cor(averaged_HOGs, method = "spearman", use = "pairwise.complete.obs")

# Convert correlation matrix to a distance matrix (1 - |correlation|)
go_dist <- as.dist(1 - go_cor_matrix)

# Perform hierarchical clustering of pathways
hc <- hclust(go_dist, method = "average")

# Cut the dendrogram to form clusters (r > 0.9 → h = 0.1)
clusters <- cutree(hc, h = 0.1) 

# Create a dataframe mapping HOGs to their cluster ID
hog_clusters <- data.frame(HOG_term = names(clusters), Cluster = clusters)

hog_t_df <- as.data.frame(t(averaged_HOGs)) # Now rows = HOG terms, cols = samples
hog_t_df$HOG_term <- rownames(hog_t_df)

# Join cluster assignments to transposed functional data
hog_t_df <- left_join(hog_t_df, hog_clusters, by = "HOG_term")

# Average the abundance of pathways within each cluster, for each sample
# (This assumes averaging makes biological sense which it does for this normalized read count data—alternatively, you could sum)
hog_cluster_avg <- hog_t_df %>%
  dplyr::select(-HOG_term) %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

# Create a label for each cluster listing all pathways it contains
hog_cluster_members <- hog_clusters %>%
  group_by(Cluster) %>%
  summarise(HOG_terms = paste(HOG_term, collapse = ", "), .groups = "drop")

# Add cluster member labels to the averaged abundance matrix
hog_cluster_avg_labeled <- as.data.frame(hog_cluster_avg) %>%
  left_join(hog_cluster_members, by = "Cluster")

# Use cluster names as rownames and remove redundant 'Cluster' column
rownames(hog_cluster_avg_labeled) <- paste0("Cluster_", hog_cluster_avg_labeled$Cluster)
hog_cluster_avg_labeled[1] <- list(NULL)

write.table(hog_cluster_avg_labeled, file = "HOG_clusters_dis.txt", sep = "\t")
write.table(hog_cluster_avg_labeled, file = "HOG_clusters_naive.txt", sep = "\t")
write.table(hog_cluster_avg_labeled, file = "HOG_clusters_AH.txt", sep = "\t")
write.table(hog_cluster_avg_labeled, file = "HOG_clusters_all.txt", sep = "\t")