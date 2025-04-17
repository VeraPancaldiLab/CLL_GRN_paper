
library(tidyverse)
library(scales)

# Define function to read and annotate enrichment files
enrichment_loader <- function(file_path) {
  file_name <- basename(file_path)
  
  # Extract metadata from file name
  time_comparison <- str_extract(file_name, "TimePoint_\\d+vs\\d+")
  patient <- str_extract(file_name, "Patient\\d+")
  database <- str_extract(file_name, "KEGG|REAC")
  gene_type <- str_extract(file_name, "DOWN|UP")
  
  # Read the enrichment file and add metadata columns
  enrichment_data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Ensure consistent data types for all columns
  enrichment_data <- enrichment_data %>%
    mutate(across(everything(), as.character))
  
  # Add metadata columns
  enrichment_data <- enrichment_data %>%
    mutate(TimeComparison = time_comparison,
           Type = gene_type,
           Patient = patient,
           Database = database)
  
  return(enrichment_data)
}

# Specify the directory containing the enrichment files
enrichment_dir_auto <- "Results/msVIPER/Autologous/Objects/"
enrichment_dir_mono <- "Results/msVIPER/Monoculture/Objects/"

# List all enrichment files
enrichment_files_auto <- list.files(enrichment_dir_auto, full.names = TRUE, pattern = "\\.csv$")
enrichment_files_mono <- list.files(enrichment_dir_mono, full.names = TRUE, pattern = "\\.csv$")

# Read and concatenate all enrichment files
enrichment_results_auto <- map_dfr(enrichment_files_auto, enrichment_loader)
enrichment_results_auto$Condition <- str_extract(enrichment_dir_auto, "Autologous|Monoculture")

enrichment_results_mono <- map_dfr(enrichment_files_mono, enrichment_loader)
enrichment_results_mono$Condition <- str_extract(enrichment_dir_mono, "Autologous|Monoculture")

enrichment_results <- rbind(enrichment_results_auto, enrichment_results_mono)

write.csv(enrichment_results, file = "Results/msVIPER/concatenated_enrichment_results_auto_mono.csv", row.names = FALSE)

time_comparison_levels <- c("TimePoint_1vs4", "TimePoint_4vs8", "TimePoint_8vs11", "TimePoint_11vs14") 

enrichment_results_filter_category <- enrichment_results %>% 
  filter(category == "Environmental Information Processing" | X == "hsa04981" | category == "Metabolism")
enrichment_results_filter_database <- enrichment_results_filter_category %>% filter(Database == "KEGG")
enrichment_results_filter_patient <- enrichment_results_filter_database %>% filter(Patient == "Patient1")

# Create dotplots for KEGG and Reactome databases
#for (db in c("KEGG", "REAC")) {
for (source in c("DOWN", "UP")) {
  # Filter data for the specific database
  db_data <- enrichment_results_filter_category %>% 
    filter(Type == source) %>%  
    mutate(TimeComparison = factor(TimeComparison, levels = time_comparison_levels))
  
  # Create dotplot with separated conditions
  dotplot_path <- ggplot(db_data, aes(x = TimeComparison, y = Description, 
                                      color = Condition, size = -log10(as.numeric(pvalue)))) +
    geom_point() +  
    facet_wrap(~Condition) +  # Separate plots for Autologous and Monoculture
    theme_minimal() +
    labs(title = paste("Enriched signalling pathways - ", source, " DA TFs"),
         x = "Time Comparison", y = "Pathway", size = "-log10(pvalue)", color = "Condition") +
    theme(axis.text.x = element_text(angle = 30, hjust = 0.5),
          strip.text = element_text(size = 16, face = "bold")) +  # Increase facet title font size and bold
    scale_y_discrete(labels = label_wrap(50)) +
    scale_color_manual(values = c("Autologous" = "#44AA99", "Monoculture" = "#E69F00"))  
  
  print(dotplot_path)  # Ensure the plot is displayed
  ggsave(paste("Results/msVIPER/enrichment_auto_mono_", source, ".svg", sep = ""), plot = dotplot_path, width = 9, height = 9)
}
#}
