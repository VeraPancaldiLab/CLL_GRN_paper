library(pheatmap)
library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(scales)
library(ggpubr)

# Load the AUCell results 
cell_of_interest <- "CLL cells (Cadot et al, 2020)"
AUCell_df <- read.csv("Supplementary_files/GRN_inference_pyScenic/Outputs/AUCELL_df.csv", header = T, row.names = 1)

# Remove the ... in the column names
colnames(AUCell_df) <- gsub("\\.\\.\\.", "", colnames(AUCell_df))

cell_names <- rownames(AUCell_df)

patient <- sub("^(P\\d)_D\\d+_.*", "\\1", cell_names)
M_values <- sub(".*_M([0-9]+)_.*", "\\1", cell_names)
D_values <- paste0("D", M_values)
day <- sub("^CPT1\\m_(M\\m+)_.*", "\\1", cell_names)

annotation_col <- data.frame(TimePoint = D_values)
rownames(annotation_col) <- cell_names

# Plot the matrix
PheatPlot <- pheatmap(t(AUCell_df),
         show_colnames = F, 
         annotation_col = annotation_col,
         fontsize = 6,
         main = paste("AUCell results: regulon activity in", cell_of_interest, " GRNBoost2 algorithm"))

# Some more filtering on the AUCell matrix
# here we seelct the top 50 (adjustable) TFs that explain the most the cell phenotypes

top_n <- 50
nonzero_mat <- AUCell_df[colSums(AUCell_df != 0) > 0, ]
top_genes <- head(order(apply(nonzero_mat, 2, var), decreasing = TRUE), top_n)
filtered_mat <- nonzero_mat[, top_genes]
TFs_ofInterest <- colnames(filtered_mat)

pheatmap::pheatmap(t(filtered_mat), 
                   scale = 'none', 
                   show_colnames = F,
                   annotation_col = annotation_col,
                   main = paste("Master TFs, ", cell_of_interest, "GRN, GRNBoost2 algorithm"),
                   fontsize = 10)#,
                   #filename = "/home/malvina.marku/Documents/Projects/RD2Bool/pySCENIC/Results/Tcell/GRNBoost/AUCell_result_topTFs_noScale.png")

