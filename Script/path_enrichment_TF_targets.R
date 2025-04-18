# Load required libraries
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(OmnipathR)
library(dplyr)
library(tidyr)
library(forcats)
library(scales)
library(ggplot2)
library(stringr)  # For text wrapping

####################################################################################################
# Load the data: 1) TF activities, 2) normalised counts
gcm_auto <- read.table("Results/Viper/Autologous/TF_score_collectri.csv", header = TRUE, row.names = 1, sep = ",")
gcm_B <- read.table("Results/Viper/B_cell/TF_score_collectri.csv", header = TRUE, row.names = 1, sep = ",")

normalized_data_auto <- read.table("Results/Normalized_data/Autologous/normalized_counts.txt", header = TRUE, row.names = 1)
normalized_data_B <- read.table("Results/Normalized_data/B_cell/normalized_counts.txt", header = TRUE, row.names = 1)

# Query Omnipath for interactions where the specified TF is the source
tf_targets <- import_all_interactions()

####################################################################################################
# Load the GRN
grn_node_features <- read.csv("Results/Temp/dynGENIE3_corr/dyngenie3_corr_full_annotation.csv", header = T)

# Define function for pathway enrichment and plotting
gprofiler_plot <- function(gene_list,
                           custom.bg,
                           gprof_sources,
                           max_set_size,
                           min_set_size,
                           correction_method,
                           n_terms_per_db,
                           title,
                           output_file) {
  
  gostres = gost(gene_list, organism = "hsapiens",
                 ordered_query = FALSE,
                 evcodes = TRUE,
                 correction_method = correction_method,
                 domain_scope = "annotated",
                 custom_bg = custom.bg,
                 sources = gprof_sources, exclude_iea = FALSE)
  if(is.null(gostres)){
    print(c("No significant enrichment results to display"))
  } else {
    gostres.df <- gostres$result
    print(paste(nrow(gostres.df), "terms / pathways significantly enriched BEFORE filtering based on minimum / maximum gene set sizes", sep =" "))
    gostres.df <- dplyr::filter(gostres$result, term_size <= max_set_size & term_size > min_set_size)
    print(paste(nrow(gostres.df), "terms / pathways significantly enriched AFTER filtering based on minimum / maximum gene set sizes", sep =" "))
    
    duplicated_levels <- duplicated(gostres.df$term_name)
    gostres.df[duplicated_levels,]$term_name <- paste(gostres.df[duplicated_levels,]$source, gostres.df[duplicated_levels,]$term_name, sep ="_")
    tmp <- as.data.frame(gostres.df %>%
                           dplyr::arrange(p_value) %>%
                           dplyr::group_by(source) %>% dplyr::slice(1:n_terms_per_db))
    tmp$term_name <- factor(tmp$term_name, levels = fct_reorder(.f = rev(tmp$term_name), .x = rev(-log10(tmp$p_value))))
    g <-  ggplot(tmp, aes(x = -log10(p_value),
                          y= term_name,
                          color = source,
                          size = recall)) +
      geom_point() +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 10)) +
      labs( x = "-log10(pvalue)",
            y = "",
            size = "Ratio",
            title = title) +
      guides(colour = guide_legend(nrow = 1, order = 1, direction = "vertical", title.position = "top", hjust =1),
             size = guide_legend(nrow = 1, order = 2, direction = "vertical", title.position = "top", hjust =-1)) + 
      scale_y_discrete(labels = label_wrap(width = 45))
    #print(g)
    # Save the plot
    ggsave(filename = output_file, plot = g, device = "svg", width = 9, height = 10)
    print(paste("Saved:", output_file))
  }
}

# Loop through Modules (1 to 11)
for (module in 1:11) {
  print(paste("Processing Module", module, "..."))
  
  # Select TFs from the current module
  tf_of_interest <- grn_node_features %>% filter(Cluster == module)
  tf_of_interest <- tf_of_interest$name
  
  # Get target genes of the TFs
  upstream_downstream_genes <- tf_targets[(tf_targets$source_genesymbol %in% tf_of_interest | 
                                             tf_targets$target_genesymbol %in% tf_of_interest), ] %>% na.omit()
  genes_of_interest <- unique(c(upstream_downstream_genes$source_genesymbol, upstream_downstream_genes$target_genesymbol))
  
  # Perform pathway enrichment & save plot
  output_path <- paste0("Results/Temp/dynGENIE3_corr/Module_enrichment/Prova_Module", module, ".svg")
  pathways_plot <- gprofiler_plot(genes_of_interest, 
                                  custom.bg = NULL, 
                                  gprof_sources = c("GO:BP", "REAC", "KEGG"), 
                                  max_set_size = 1000, 
                                  min_set_size = 10, 
                                  correction_method = 'fdr', 
                                  n_terms_per_db = 10, 
                                  title = paste("Pathway Enrichment, Module", module),
                                  output_file = output_path)
}

print("All modules processed successfully.")