## ---------------------------------------------------------------------------------------------------------------------------------------
library(reshape2)
library(doRNG)
library(doParallel)
library(sva)
library(tidyverse)
library(DESeq2)
library(clusterSim)
library(dplyr)
library(ggplot2)
library(factoextra) 
library(FactoMineR)
library(tibble)
library(qpcR)
library(MASS)
library(gridExtra)


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: remove_ensembl_id
# Description:
#   Removes the ENSEMBL Gene ID prefix from the row names of a gene count matrix, keeping only the Gene Symbol.
#   If duplicate Gene Symbols are found, it appends a numeric suffix to differentiate them.
#
# Arguments:
#   gene_count_matrix - A matrix or data frame where row names include ENSEMBL IDs followed by an underscore and Gene Symbols.
#
# Returns:
#   The input gene count matrix with modified row names containing only Gene Symbols (with suffixes for duplicates if necessary).
remove_ensembl_id <- function(gene_count_matrix) {
  rownames_original <- rownames(gene_count_matrix)
  new_rownames <- sub(".*_", "", rownames_original)
  
  if (any(duplicated(new_rownames))) {
    suffix <- ave(new_rownames, new_rownames, FUN = seq_along)
    new_rownames[duplicated(new_rownames)] <- paste0(new_rownames[duplicated(new_rownames)], "_", suffix[duplicated(new_rownames)])
  }
  
  rownames(gene_count_matrix) <- new_rownames
  return(gene_count_matrix)
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: create_patient_replicate_tables
# Description:
#   Creates separate tables for each patient and replicate combination by selecting columns
#   that contain a specific pattern in their names from a given gene count matrix.
#
# Arguments:
#   gene_count_matrix - A data frame or matrix of gene counts.
#   patient_count     - A vector indicating the patient identifiers (e.g., c(1, 2, 3)).
#   replicate_count   - A vector indicating the replicate identifiers (e.g., c(1, 2)).
#   matrix_type       - A character string that defines the type of matrix ("auto" or "B").
#
# Returns:
#   This function does not return a value directly; it creates a list of patient-specific tables,
#   assigns them to a list 'patient_tables', and also assigns each table to the global environment.
create_patient_replicate_tables <- function(gene_count_matrix, patient_count, replicate_count, matrix_type) {
  patient_tables <- list()
  matrix_type_str <- ""
  
  if (identical(matrix_type, "hetero")) {
    matrix_type_str <- "_hetero"
  } else if (identical(matrix_type, "auto")) {
    matrix_type_str <- "_auto"
  } else if (identical(matrix_type, "B")) {
    matrix_type_str <- "_B"
  } else {
    stop("Invalid matrix_type provided. Please use a matrix corresponding to a valid condition.")
  }
  
  for (i in patient_count) {
    for (j in replicate_count) {
      table_name <- paste0("Patient_", i, "_rep_", j, matrix_type_str, "_table")
      patient_table <- gene_count_matrix %>% dplyr::select(contains(paste0("_p", i, "_rep", j)))
      patient_tables[[table_name]] <- patient_table
      assign(table_name, patient_table, envir = .GlobalEnv)
    }
  }
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: apply_impute_mean
# Description:
#   Applies a mean imputation function (impute.mean) to each column of patient-specific gene count tables,
#   updating the tables in the global environment and printing the number of infinite values per table.
#
# Arguments:
#   patient_count - A vector of patient identifiers.
#   replicate_count - A vector of replicate identifiers.
#   conditions - A vector of conditions (e.g., "auto", "B") that are part of the table name.
#
# Returns:
#   No explicit return value. The function updates tables in the global environment and prints status messages.
apply_impute_mean <- function(patient_count, replicate_count, conditions) {
  for (i in patient_count) {
    for (j in replicate_count) {
      for (condition in conditions) {
        table_name <- paste0("Patient_", i, "_rep_", j, "_", condition, "_table")
        patient_table <- get(table_name, envir = .GlobalEnv)
        
        # Apply impute.mean to each column and update the patient_table
        patient_table <- apply(patient_table, 2, impute.mean)
        
        # Assign the updated table back to the global environment
        assign(table_name, patient_table, envir = .GlobalEnv)
        
        # Print the number of infinite values in the updated table
        inf_values <- sum(apply(patient_table, 2, function(x) sum(is.infinite(x))))
        cat(paste0("Number of infinite values in ", table_name, ": ", inf_values, "\n"))
      }
    }
  }
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: process_patient_tables
# Description:
#   Processes patient replicate tables for each condition by filtering rows with zero variance,
#   identifying common genes across all tables, and creating a combined time series data matrix.
#
# Arguments:
#   patient_count - A vector of patient identifiers.
#   replicate_count - A vector of replicate identifiers.
#   conditions - A vector of conditions to process (e.g., "auto", "B").
#
# Returns:
#   A list containing:
#     common_genes    - A list of common gene names for each condition and across all conditions.
#     time_series_data - A combined matrix of filtered patient tables with common genes.
process_patient_tables <- function(patient_count, replicate_count, conditions) {
  common_genes <- list()
  time_series_data <- list()
  
  # Iterate through conditions and process tables
  for (condition in conditions) {
    # Initialize empty lists to store filtered tables
    filtered_tables <- list()
    
    for (i in patient_count) {
      for (j in replicate_count) {
        table_name <- paste0("Patient_", i, "_rep_", j, "_", condition, "_table")
        patient_table <- get(table_name, envir = .GlobalEnv)
        
        # Add a row indicating the time point
        #time_points <- c(1:ncol(patient_table))
        #patient_table <- rbind(time_points, patient_table)
        
        # Remove rows with 0 variance
        patient_table <- patient_table[rowSums(patient_table[]) > 0,]
        
        # Store filtered table
        filtered_tables[[paste0("Patient_", i, "_rep_", j)]] <- patient_table
      }
    }
    
    # Find common genes across all filtered tables
    common_genes[[condition]] <- Reduce(intersect, lapply(filtered_tables, rownames))
    
    # Keep only the genes that appear in common_genes
    for (table_name in names(filtered_tables)) {
      filtered_table <- filtered_tables[[table_name]]
      filtered_table <- as.matrix(subset(filtered_table, rownames(filtered_table) %in% common_genes[[condition]]))
      
      # Assign the updated table back to the global environment
      assign(paste0(table_name, "_", condition, "_filter"), filtered_table, envir = .GlobalEnv)
      
      # Add filtered_table to time_series_data list
      time_series_data[[paste0(table_name, "_", condition)]] <- filtered_table[1:nrow(filtered_table),]
    }
  }
  
  # Combine filtered tables to create the final time_series_data list
  time_series_data <- do.call(cbind, time_series_data)
  colnames(time_series_data) <- gsub(".*\\.", "", colnames(time_series_data))
  
  # Combine common genes from all conditions
  all_common_genes <- Reduce(intersect, common_genes)
  
  return(list(common_genes = all_common_genes, time_series_data = time_series_data))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: process_time_series_data
# Description:
#   Processes time series gene count data by creating a DESeq2 dataset,
#   filtering low count genes, normalizing counts, and removing genes with zero variance.
#
# Arguments:
#   time_series_data - A matrix of gene counts where columns represent samples.
#
# Returns:
#   A list containing:
#     norm   - A data frame of normalized counts with non-zero variance.
#     dds_obj - The DESeq2 dataset object after filtering and size factor estimation.
#     unnorm - The unnormalized counts adjusted by the estimated size factors.
process_time_series_data <- function(time_series_data) {
  
  sample_info = as.data.frame(stringr::str_split_fixed(colnames(time_series_data), pattern = "_", n = 3))
  names(sample_info)[names(sample_info) == "V1"] <- "Day"
  names(sample_info)[names(sample_info) == "V2"] <- "Patient"
  names(sample_info)[names(sample_info) == "V3"] <- "Replicate"
  
  # Create DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(time_series_data), 
                                        colData = sample_info, 
                                        design = ~ Day + Patient + Replicate)
  #nrow(dds)
  
  # Remove rows having low counts
  firstcondition <- rowSums(counts(dds) > 1) >= 2 # at least 2 samples with a count of 2 or higher
  dds <- dds[firstcondition,]
  #nrow(dds)
  
  # Generate normalized counts
  dds <- DESeq2::estimateSizeFactors(dds)
  sizeFactors(dds)
  normalized_counts <- counts(dds, normalized = TRUE)
  unnormalized_counts <- t(t(normalized_counts) * sizeFactors(dds))
  
  # Remove all genes with zero variance (i.e constant expression)
  variances = matrixStats::rowVars(normalized_counts)
  range(variances)
  normalized_counts = normalized_counts[!matrixStats::rowVars(normalized_counts) == 0,] %>% as.data.frame()
  
  return(list(norm = normalized_counts, dds_obj = dds, unnorm = unnormalized_counts))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: plot_PCA
# Description:
#   Performs a Principal Component Analysis (PCA) on the provided data and plots the first two principal components.
#   The function extracts patient information from the data's column names for labeling.
#
# Arguments:
#   data      - A numeric matrix or data frame where columns are samples.
#   data_type - A descriptive label for the type of data (e.g., "Normalized counts").
#   condition - A label describing the experimental condition.
#
# Returns:
#   A ggplot object representing the PCA plot. 
plot_PCA <- function(data, data_type, condition) {
  
  # Address high explained variance in PC1
  #data_log <- log1p(data)
  
  # Perform PCA
  PCA <- prcomp(t(data), scale = FALSE)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # Calculate ratio of standard deviations for axis scaling
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  # Create df with PCA coordinates and sample names
  dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
  
  # Extract patient information from column names
  dataGG <- dataGG %>% rownames_to_column("Patient") %>% 
    mutate(P = str_split_fixed(colnames(data), pattern = "_", n = 3)[, 2],
           J = str_split_fixed(colnames(data), pattern = "_", n = 3)[, 1])
  
  # Plot PCA with labeled points
  title <- paste0("PCA - ", data_type, " - ", condition)
  ggplot(dataGG, aes(PC1, PC2, label = Patient, col = P)) +
    geom_point() +
    geom_text(aes(label = Patient), hjust = -.2, vjust = .4) +
    ggtitle(title) + 
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: patient_pca
# Description:
#   Generates PCA plots for each patient separately using DESeq2 sample information and normalized counts.
#   It uses FactoMineR for PCA computation and factoextra for visualization.
#
# Arguments:
#   dds                - A DESeq2 dataset object containing sample information.
#   normalized_counts  - A matrix or data frame of normalized counts.
#   condition          - A character string describing the experimental condition (e.g., "Autologous" or "B cell").
#
# Returns:
#   No explicit return value. The function displays PCA plots for each patient arranged in a grid.
patient_pca <- function(dds, normalized_counts, condition) {
  
  unique_patients <- unique(dds$Patient)
  patient_plots <- list()
  
  for (patient in unique_patients) {
    res.pca <- FactoMineR::PCA(
      cbind(
        data.frame("condition" = dds$Replicate, "Patient" = dds$Patient, "Day" = dds$Day),
        t(normalized_counts)
      ) %>% dplyr::filter(Patient == patient),
      graph = FALSE,
      quali.sup = c(1:3)
    )
    
    patient_plot <- factoextra::fviz_pca_ind(
      res.pca,
      habillage = "Day",
      repel = TRUE,
      pointsize = 4,
      mean.point = FALSE
    ) + ggtitle(paste0("Individuals PCA - ", condition, " - Patient ", patient))
    
    patient_plots[[paste0("Patient_", patient)]] <- patient_plot
  }
  
  gridExtra::grid.arrange(grobs = patient_plots, nrow = 1)
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: process_normalized_counts
# Description:
#   Processes normalized counts to extract only the rows corresponding to TFs.
#   The function uses a provided TF names file and filters the counts matrix based on the chosen identifier form.
#
# Arguments:
#   normalized_counts - A data frame or matrix of normalized gene counts.
#   tf_names_file     - A file path to a CSV file containing human TF names.
#   id_form           - A character string specifying the identifier format ("genenames" or "ensemblIDs_genenames").
#   sample_suffixes   - A vector of suffix strings to subset the normalized counts for different samples.
#
# Returns:
#   A list containing:
#     normalized_counts_TF - A data frame of normalized counts for the selected TFs.
#     normalized_counts_list - A list of matrices corresponding to different patients/replicates.
#     TFs_vector           - A vector of TF names used for filtering.
process_normalized_counts <- function(normalized_counts, tf_names_file, id_form, sample_suffixes) {
  
  # Keep only rows corresponding to TF
  normalized_counts <- as.data.frame(normalized_counts)
  
  if (id_form == "genenames") {
    TFs_df <- read.csv(tf_names_file, header = FALSE)
    TFs_vector <- TFs_df[['V1']]
    TFs_vector <- TFs_vector[TFs_vector %in% rownames(normalized_counts)] 
    normalized_counts_TF <- normalized_counts[rownames(normalized_counts) %in% TFs_vector, ]
  } else if (id_form == "ensemblIDs_genenames") {
    TFs_df <- read.csv(tf_names_file, header = FALSE)
    TFs_vector <- TFs_df[['V1']]
    grepl(paste(TFs_vector, collapse = "|"), rownames(normalized_counts))
    TFs_vector <- TFs_vector[TFs_vector %in% (str_extract(rownames(normalized_counts), paste(TFs_vector, collapse = "|")))]
    normalized_counts_TF <- normalized_counts[grepl(paste0(TFs_vector, collapse = "|"), rownames(normalized_counts)), ]
  }
  
  # List of matrices corresponding to the different patients/replicates 
  df_conditions <- sapply(sample_suffixes,
                          function(x) normalized_counts_TF[endsWith(names(normalized_counts_TF), x)],
                          simplify = FALSE)
  
  normalized_counts_list <- lapply(df_conditions, function(x) as.matrix(x))
  
  # Return the results as a list
  return(list(normalized_counts_TF = normalized_counts_TF, normalized_counts_list = normalized_counts_list, TFs_vector = TFs_vector))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: sample_pca
# Description:
#   Constructs a PCA plot for samples based on input data.
#   Extracts sample metadata from column names and visualizes the first two principal components.
#
# Arguments:
#   input_data - A numeric matrix or data frame of gene expression data where columns represent samples.
#   Culture    - A character string describing the culture type (e.g., "Autologous", "B cell").
#
# Returns:
#   A ggplot object representing the PCA plot for the samples.
sample_pca <- function(input_data, Culture){
  
  pca.res <- prcomp(t(input_data)) #, scale. = TRUE
  scores = as.data.frame(pca.res$x)
  
  percentVar <- round(100 * pca.res$sdev ^ 2 / sum(pca.res$sdev ^ 2), 1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  dataGG <- data.frame(PC1 = pca.res$x[, 1], PC2 = pca.res$x[, 2])
  dataGG <- dataGG %>% rownames_to_column("Patient") %>% mutate(
    Day = str_split_fixed(colnames(input_data), pattern = "_", n = 4)[, 1] %>% as.factor(),
    Patient = str_split_fixed(colnames(input_data), pattern = "_", n = 4)[, 2],
    Replicate = str_split_fixed(colnames(input_data), pattern = "_", n = 4)[, 3],
    Condition = str_split_fixed(colnames(input_data), pattern = "_", n = 4)[, 4])
  dataGG$Day <- factor(dataGG$Day, levels = c("D1", "D4", "D8", "D11", "D14"))
  
  ggplot(dataGG, aes(PC1, PC2, label = Day, col = Condition, shape = Patient)
  ) +
    geom_point() +
    geom_text(aes(label = Day), hjust = -.25, vjust = -.25) +
    ggtitle(paste("PCA plot:", Culture)) +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    expand_limits(x = 100, y = 100) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    scale_color_manual(values = c("auto" = "#44AA99", "B" = "#E69F00"))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: plot_day_data
# Description:
#   Visualizes data for each day and condition by melting the data into a long format,
#   calculating averages and standard deviations, and plotting both individual patient trends and overall means.
#
# Arguments:
#   data - A data frame containing day information and measurements for different patients and conditions.
#
# Returns:
#   No explicit return value. The function generates and prints a ggplot showing trends over days.
plot_day_data <- function(data) {
  
  # Extract day info from the first column (Day)
  data$Day <- as.numeric(str_extract(data$Day, "(?<=D)[0-9]+"))
  
  # Melt the dataframe to long format for easier plotting
  data_long <- data %>%
    pivot_longer(
      cols = starts_with("P"),
      names_to = c("Patient", "Condition"),
      names_pattern = "P(\\d+)_(\\w+)",
      values_to = "Value"
    )
  
  data_long_p1 <- data_long %>% filter(Patient == 1)
  data_long_p2 <- data_long %>% filter(Patient == 2)
  data_long_p3 <- data_long %>% filter(Patient == 3)
  
  data_long <- data_long %>%
    mutate(Day = as.numeric(Day),  # Use numeric Day for continuous lines
           Patient = as.factor(Patient),
           Condition = as.factor(Condition),
           Value = as.numeric(Value))
  
  # Calculate the average value for each patient across replicates
  patient_avg <- data_long %>%
    group_by(Day, Patient, Condition) %>%
    summarise(Avg_Value = mean(Value, na.rm = TRUE)) %>%
    ungroup()
  
  # Calculate means and standard error for each condition and day
  summary_data <- data_long %>%
    group_by(Day, Condition) %>%
    summarise(Mean = mean(Value, na.rm = TRUE),
              SD = sd(Value, na.rm = TRUE)) %>%
    ungroup()
  
  # Plot the values
  p <- ggplot() +
    # Dashed line for the average of each patient across replicates
    geom_line(data = patient_avg, aes(x = Day, y = Avg_Value, color = Condition, group = interaction(Patient, Condition)), 
              linetype = "dashed", size = 0.7) +
    # Solid line for the mean values across all patients
    geom_line(data = summary_data, aes(x = Day, y = Mean, color = Condition, group = Condition), size = 1.2) +
    # Error bars for the mean values
    geom_errorbar(data = summary_data, aes(x = Day, y = Mean, ymin = Mean - SD, ymax = Mean + SD, color = Condition), 
                  width = 0.2, size = 0.5) +
    geom_point(data = data_long, aes(x = Day, y = Value, color = Condition, shape = Patient), 
               position = position_jitter(width = 0.2), size = 3) +
    labs(title = "CLL cell viability", x = "Days", y = "Viability (%)") +
    scale_color_manual(values = c("auto" = "#44AA99", "B" = "#E69F00")) +
    theme_minimal() +
    theme(legend.position = "top")
  
  print(p)
}

