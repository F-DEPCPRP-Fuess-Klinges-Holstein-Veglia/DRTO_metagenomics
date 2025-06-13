library(tidyverse)

func_cluster = read.table(file = "functional_pathway_clusters_all.txt", header = TRUE, )
hog_cluster = read.table(file = "HOG_clusters_all.txt", header = TRUE, )
#meta = read.csv(file = "metadata.csv", header = TRUE, )

func_cluster <- func_cluster[ , grep("\\.Dis$", names(func_cluster))]
hog_cluster <- hog_cluster[ , grep("\\.Dis$", names(hog_cluster))]

#func_cluster <- func_cluster[ , grep("\\.N$", names(func_cluster))]
#hog_cluster <- hog_cluster[ , grep("\\.N$", names(hog_cluster))]

rownames(hog_cluster) <- paste0("OG_", rownames(hog_cluster))
rownames(func_cluster) <- paste0("Pathway_", rownames(func_cluster))

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

sig_results <- results_long %>%
  filter(abs(p_adj) < 0.05) %>%
  filter(abs(Correlation) > 0.79) %>%
  # only keep significant correlations
  arrange(desc(abs(Correlation)))
length(sig_results$p_adj)
length(unique(sig_results$HOG_ID))
length(unique(sig_results$Pathway))

write.table(filtered_results, file = "high_corr_AH.txt", sep= "\t")
write.table(filtered_results, file = "high_corr_disease.txt", sep= "\t")
write.table(filtered_results, file = "high_corr_naive.txt", sep= "\t")
write.table(sig_results, file = "high_corr_dis_all.txt", sep= "\t")

write.table(sig_results, file = "sig_corr_dis.txt", sep= "\t")
