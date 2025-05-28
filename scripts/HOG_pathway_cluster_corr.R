library(tidyverse)

func_cluster = read.table(file = "functional_pathway_clusters_dis.txt", header = TRUE, )
hog_cluster = read.table(file = "HOG_clusters_dis.txt", header = TRUE, )
#meta = read.csv(file = "metadata.csv", header = TRUE, )

rownames(hog_cluster) <- paste0("OG_", rownames(hog_cluster))
rownames(func_cluster) <- paste0("Pathway_", rownames(func_cluster))

#exporting metadata as a separate file to make dataframe purely numeric
func_meta <- as.data.frame(func_cluster$Pathway)
rownames(func_meta) <- rownames(func_cluster)
#write.table(func_meta, file = "func_clusters_names_dis.txt", sep = "\t")
func_cluster <- func_cluster %>% select(-Pathway)
rm(func_meta)

hog_meta <- as.data.frame(hog_cluster$HOG_terms)
rownames(hog_meta) <- rownames(hog_cluster)
#write.table(hog_meta, file = "hog_clusters_names.txt", sep = "\t")
hog_cluster <- hog_cluster %>% select(-HOG_terms)
rm(hog_meta)

func <- as.data.frame(t(func_cluster))
rm(func_cluster)
hogs <- as.data.frame(t(hog_cluster))
rm(hog_cluster)

n_hog <- ncol(hogs)
n_pathway <- ncol(func)

cor_mat <- matrix(NA, nrow = n_hog, ncol = n_pathway)
pval_mat <- matrix(NA, nrow = n_hog, ncol = n_pathway)

# Assign row/column names
rownames(cor_mat) <- colnames(hogs)
colnames(cor_mat) <- colnames(func)
rownames(pval_mat) <- colnames(hogs)
colnames(pval_mat) <- colnames(func)

# Loop over all GO x pathway combinations
for (i in seq_len(n_hog)) {
  for (j in seq_len(n_pathway)) {
    x <- hogs[[i]]
    y <- func[[j]]
    
    # Use cor.test to get both correlation and p-value
    result <- cor.test(x, y, method = "spearman", exact = FALSE)
    
    cor_mat[i, j] <- result$estimate
    pval_mat[i, j] <- result$p.value
  }
}

results_long <- as_tibble(cor_mat, rownames = "HOG_ID") %>%
  pivot_longer(-HOG_ID, names_to = "Pathway", values_to = "Correlation") %>%
  mutate(p_value = as.vector(pval_mat)) %>%
  arrange(p_value)
results_long <- results_long %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

filtered_results <- results_long %>%
  filter(abs(Correlation) > 0.79) %>%  # keep strong positive and negative correlations
  arrange(desc(abs(Correlation))) 

sig_results <- results_long %>%
  filter(abs(p_adj) < 0.05) %>%  # only keep significant correlations
  arrange(desc(abs(Correlation)))
length(sig_results$p_adj)
length(unique(sig_results$HOG_ID))
length(unique(sig_results$Pathway))

write.table(filtered_results, file = "high_corr_AH.txt", sep= "\t")
write.table(filtered_results, file = "high_corr_disease.txt", sep= "\t")
write.table(filtered_results, file = "high_corr_naive.txt", sep= "\t")