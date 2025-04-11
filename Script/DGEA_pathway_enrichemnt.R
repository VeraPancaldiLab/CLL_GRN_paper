library(DESeq2)
library(gprofiler2)
library(dplyr)
library(forcats)
library(ggplot2)
library(tidyverse)
library(scales)
library(ggpubr)

####################################################################################################
# Load the data: non-normalized counts

data_auto <- read.table("Results/Normalized_data/Autologous/unnormalized_counts.txt", 
                                   header = TRUE, row.names = 1) 
colnames(data_auto) <- gsub("^_auto", "autologous", colnames(data_auto))
data_B <- read.table("Results/Normalized_data/B_cell/unnormalized_counts.txt", 
                                header = TRUE, row.names = 1) 
colnames(data_B) <- gsub("^_B", "monoculture", colnames(data_B))

####################################################################################################
# Combine the two conditions into a single dataframe
# Select the common genes between the two datasets and merge them

common_genes <- intersect(rownames(data_auto), rownames(data_B))

data_auto_filter <- data_auto[common_genes, , drop = FALSE] %>% drop_na()
data_B_filter <- data_B[common_genes, , drop = FALSE] %>% drop_na()

combined_data <- cbind(data_auto_filter, data_B_filter)

####################################################################################################
# Run this part if we want to do the analysis for each patient separately

combined_data_P1 <- combined_data %>% select(contains("_p1_"))
combined_data_P2 <- combined_data %>% select(contains("_p2_"))
combined_data_P3 <- combined_data %>% select(contains("_p3_"))

####################################################################################################
# Get the information from the combined dataset

combined_data_J1 <- data_auto %>% dplyr::select(starts_with("J1_"))
combined_data_J4 <- data_auto %>% dplyr::select(starts_with("J4_"))
combined_data_J8 <- data_auto %>% dplyr::select(starts_with("J8_"))
combined_data_J11 <- data_auto %>% dplyr::select(starts_with("J11_"))
combined_data_J14 <- data_auto %>% dplyr::select(starts_with("J14_"))

####################################################################################################

sample_info = as.data.frame(stringr::str_split_fixed(colnames(data_B), pattern = "_", n = 4))
names(sample_info)[names(sample_info) == "V1"] <- "Day"
names(sample_info)[names(sample_info) == "V2"] <- "Patient"
names(sample_info)[names(sample_info) == "V3"] <- "Replicate"
names(sample_info)[names(sample_info) == "V4"] <- "Condition"

ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = round(data_B),
                                            colData = sample_info,
                                            design = ~ Day)
dds = DESeq(ddsMat)
          # Run this lines only if there are not enough samples for statistics 
          dds <- estimateSizeFactors(ddsMat)
          dds <- estimateDispersionsGeneEst(dds)
          dispersions(dds) <- mcols(dds)$dispGeneEst
          dds <- nbinomWaldTest(dds)

contrast_condition <- c("Day", "J4", "J1")
results = results(dds, contrast = contrast_condition)

# keep only the significant genes
results_sig = subset(results, padj < 0.05)
# get the significant up-regulated genes
up = subset(results_sig, log2FoldChange > 1.5)
up_df <- as.data.frame(up)
up_df <- up_df[order(up_df$padj), ]
write.csv(up_df, "Results/Differential_expression/Auto_mono/Mono_D4_vs_D1_up.csv", quote = F)
# get the significant down-regulated genes
down = subset(results_sig, log2FoldChange < -1.5)
down_df <- as.data.frame(down)
down_df <- down_df[order(down_df$padj), ]
write.csv(down_df, "Results/Differential_expression/Auto_mono/Mono_D4_vs_D1_down.csv", quote = F)

##################################################################################

# Perform enrichment
gprofiler_plot <- function(gene_list,
                           custom.bg,
                           gprof_sources,
                           max_set_size,
                           min_set_size,
                           correction_method,
                           n_terms_per_db,
                           title){
  
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
    print(g)
    # png(filename = paste(Figures_dir, gsub(pattern = ".xlsx", replacement = ".png", output_file_name), sep ="/"), width = 300, height = 300, units = "cm", res = 300)
    # g
    # dev.off()
  }
}

# Use the up- or down- genes for the enrichment
pathways_plot <- gprofiler_plot(up@rownames, 
                        custom.bg = NULL, 
                        gprof_sources = c("GO:BP", "REAC", "KEGG"), 
                        max_set_size = 1000, 
                        min_set_size = 10, 
                        correction_method = 'fdr', n_terms_per_db = 10, title = "Monoculture D4 vs D1, up")
ggsave("Results/Differential_expression/Auto_mono/pathways_plot_mono_D4_vs_D1_up.svg", plot = pathways_plot, width = 10, height = 9, dpi = 1000)
####################################################################################################
# volcano plot

results$significance <- ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 1.5, "Significant", "Not Significant")
results$Gene <- rownames(results)  # Add gene names for labeling

# Convert results to a data frame for plotting
res_df <- as.data.frame(results)
# Filter for significant genes (adjusted p-value < 0.05 and |log2FoldChange| > 1)
significant_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1.5, ]

# Rank by adjusted p-value and select the top 15 significant genes
top_genes <- head(significant_genes[order(significant_genes$padj), ], 20)

# Add a column to mark these top 15 genes for labeling
res_df$top_label <- ifelse(res_df$Gene %in% top_genes$Gene, res_df$Gene, NA)

# Count the number of significant up- and down-regulated genes
upregulated_count <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1.5, na.rm = TRUE)
downregulated_count <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1.5, na.rm = TRUE)


# Create the volcano plot with gene names for the top 15 significant genes
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +  # Plot points
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "#7B5EA3")) +
  theme_minimal() +
  
  # Add dashed line for significance threshold (p-value = 0.05, log fold change threshold = ±1)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # p-value threshold line
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "blue") +
  
  
  labs(title = "DGEA: autologous vs monoculture",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.title = element_blank()) +
  
  geom_text(aes(label = top_label), vjust = 1.5, hjust = 0.5, size = 3, color = "black", na.rm = TRUE) +
  
  # Add text labels for the top 15 genes
  geom_text(aes(label = top_label), vjust = 1.5, hjust = 0.5, size = 3, color = "black", na.rm = TRUE)

######################################################################################################

# Build a MA plot (log ratio - mean average)
results_df <- as.data.frame(results)
ggma_plot <- ggmaplot(results_df, main = "Monoculture, D4 vs D1",
         fdr = 0.05, fc = 1.5, 
         size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(results_df$name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())
ggsave("Results/Differential_expression/Auto_mono/MA_plot_mono_D4_vs_D1.svg", plot = ggma_plot, width = 6, height = 5, dpi = 600)
