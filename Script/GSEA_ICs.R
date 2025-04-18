# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Change if using another organism
library(DOSE)
library(ggplot2)
library(dplyr)
library(forcats)
library(ReactomePA)
library(stringr)  # For wrapping text

# Load gene list with scores
ICscore <- read.table("Results/ICA/work/combined_auto_mono_after_normalized_df_ICA/combined_auto_mono_after_normalized_df_ica_S.xls", header = TRUE)

# Loop through IC1 to IC10
for (i in 1:10) {
  IC_name <- paste0("IC", i)
  print(paste("Processing", IC_name, "..."))  # Progress message
  
  # Create a gene list for the current IC
  gene_list <- data.frame(PROBE = ICscore$PROBE, Score = ICscore[[IC_name]])
  gene_list <- gene_list[order(gene_list$Score, decreasing = TRUE), ]
  
  # Convert gene symbols to ENTREZ IDs
  gene_list$entrez_id <- mapIds(org.Hs.eg.db, 
                                keys = as.character(gene_list$PROBE), 
                                column = "ENTREZID", 
                                keytype = "SYMBOL", 
                                multiVals = "first")
  
  # Remove NA values
  gene_list <- gene_list[!is.na(gene_list$entrez_id), ]
  
  # Create a named numeric vector for GSEA
  gene_vector <- setNames(gene_list$Score, gene_list$entrez_id)
  #all_genes <- unique(mapIds(org.Hs.eg.db, keys = ICscore$PROBE, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  
  # Perform GSEA for REACTOME pathways
  gsea_result <- gsePathway(geneList = gene_vector,
                            organism = "human",
                            minGSSize = 10,
                            maxGSSize = 2000,
                            pvalueCutoff = 0.05)
  
  # Perform GSEA for KEGG pathways
  # gsea_result_kegg <- gseKEGG(geneList = gene_vector,
  #                             organism = "hsa",
  #                             keyType = "ncbi-geneid",  # Correct KEGG key type
  #                             minGSSize = 10,
  #                             maxGSSize = 1000,
  #                             pvalueCutoff = 0.05)
  
  # Perform GSEA for GO pathways
  # gsea_result_go <- gseGO(geneList = gene_vector,
  #                         OrgDb = org.Hs.eg.db,
  #                         keyType = "ENTREZID",
  #                         minGSSize = 10,
  #                         maxGSSize = 1000,
  #                         pvalueCutoff = 0.05)
  
  # Convert results to dataframe
  gsea_df <- as.data.frame(gsea_result)  # Choose database results for visualization
  
  # Filter significant pathways
  gsea_filtered <- gsea_df %>% filter(p.adjust <= 0.05)
  
  # Check if there are significant pathways before proceeding
  if (nrow(gsea_filtered) > 0) {
    # Select top 10 positively and negatively enriched pathways
    top_up <- gsea_filtered %>% arrange(desc(NES)) %>% head(17)
    top_down <- gsea_filtered %>% arrange(NES) %>% head(17)
    
    # Combine selected pathways
    top_terms <- bind_rows(top_up, top_down)
    
    # Define color based on NES sign
    top_terms$NES_category <- ifelse(top_terms$NES > 0, "Positive NES", "Negative NES")
    
    # Wrap pathway names for better readability
    top_terms$Description <- str_wrap(top_terms$Description, width = 40)  # Adjust width if needed
    
    # Create the plot
    p <- ggplot(top_terms, aes(x = NES, y = fct_reorder(Description, NES), fill = NES_category)) +
      geom_col() +
      scale_fill_manual(values = c("Positive NES" = "#4c7fa5", "Negative NES" = "#a2893c"), guide = "none") +
      labs(x = "Normalized Enrichment Score", y = NULL, title = paste("GSEA Reactome Enrichment,", IC_name)) +
      theme_minimal() +
      theme(legend.title = element_blank(), text = element_text(size = 16))
    
    # Save the plot as a high-resolution figure
    filename <- paste0("../CLL_GRN_paper/Results/ICA/GSEA_ICs/GSEA_Reactome_", IC_name, ".png")
    ggsave(filename, plot = p, width = 10, height = 15, dpi = 600, bg = "white")  # High-resolution
  } else {
    print(paste("No significant pathways found for", IC_name, "(p.adjust ≤ 0.05)."))
  }
}

print("All ICs processed successfully.")
