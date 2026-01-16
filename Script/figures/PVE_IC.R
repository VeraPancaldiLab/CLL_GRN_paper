normalized_counts_auto <- read.csv("Results/Normalized_data/Autologous/autologous_log2normalized_counts.txt", header = T, row.names = 1)

normalized_counts_B <- read.csv("Results/Normalized_data/B_cell/monoculture_log2normalized_counts.txt", header = T, row.names = 1)


# Filter the count matrices according to the gene list of interest and create a common dataframe with the two conditions
common_genes_auto_mono <- intersect(rownames(normalized_counts_auto), rownames(normalized_counts_B))

normalized_counts_auto <- normalized_counts_auto[common_genes_auto_mono,  , drop = F]
normalized_counts_B <- normalized_counts_B[common_genes_auto_mono,  , drop = F]

combined_df <- cbind(normalized_counts_auto, normalized_counts_B)

metagene_matrix <- read.csv("Results/ICA/work/combined_auto_mono_after_normalized_df_ICA/combined_auto_mono_after_normalized_df_ica_S.xls", 
                            header = T, 
                            row.names = 1,
                            sep = "\t")
metagene_matrix$X <- NULL

common_genes_ICA_data <- intersect(rownames(combined_df), rownames(metagene_matrix))

data_cut <- combined_df[common_genes_ICA_data, , drop = FALSE]

compute_variance_explained <- function(expression_matrix, metagene_matrix, IC_index) {
  # Extract the chosen Independent Component's loadings (gene scores)
  IC <- metagene_matrix[, IC_index]
  
  # Project gene expression onto the IC using linear regression
  fitted_values <- tcrossprod(IC, qr.solve(IC, expression_matrix))
  
  # Compute total variance in the original gene expression matrix
  total_variance <- sum(apply(expression_matrix, 1, var))
  
  # Compute variance explained by the IC
  explained_variance <- sum(apply(fitted_values, 1, var))
  
  # Compute proportion of variance explained
  PVE <- 100*explained_variance / total_variance
  
  return(PVE)
}

num_ICs <- ncol(metagene_matrix)

# Compute PVE for all ICs
PVE_values <- sapply(1:num_ICs, function(IC_index) {
  compute_variance_explained(data_cut, metagene_matrix, IC_index)
})

# Plot results
barplot(PVE_values, names.arg = paste0("IC", 1:num_ICs), 
        main = "PVE by each IC",
        xlab = "PVE", ylab = "Independent Component",
        col = "steelblue", las = 1, horiz = T)
png(filename = paste0("Results/ICA/Figures/PVE.png"), width = 3, height = 10)
ggsave(filename = paste0("Results/ICA/Figures/PVE.svg"), plot = p, device = "svg", width = 3, height = 10, , dpi = 600)

################################################################################
# Get the list of the genes with the most extreme loadings in a given IC to investigate further 

top_genes <- rownames(metagene_matrix %>% filter(metagene_matrix$IC3 > 2))

pathways_plot <- gprofiler_plot(top_genes, 
                                custom.bg = NULL, 
                                gprof_sources = c("GO:BP", "REAC", "KEGG"), 
                                max_set_size = 1000, 
                                min_set_size = 10, 
                                correction_method = 'fdr', n_terms_per_db = 10, title = "IC1, positive")
