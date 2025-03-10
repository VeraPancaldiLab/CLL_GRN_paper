# Biodica.R
# This file contains all functions used in the BIODICA analysis pipeline.
# Each function includes a description, its arguments, and its return value.

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(gridExtra)
library(lsa)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dorothea)
library(viper)
library(ggplot2)
library(dendextend)
library(data.table)

# ------------------------------------------------------------------------------
# Function: merge_replicates
# Description: Merges replicate rows in a combined metasample data frame by averaging 
#              the replicates for each unique sample (based on a sample prefix).
#
# Arguments:
#   metasample_df_combined: A data frame containing the combined metasample data
#                           with replicates in the row names.
#
# Returns:
#   A data frame where replicates have been merged (averaged) for each unique sample.
# ------------------------------------------------------------------------------
merge_replicates <- function(metasample_df_combined) {
  # Extract unique sample prefixes by removing replicate labels from row names.
  sample_prefixes <- unique(gsub("_rep[12]_(auto|B)", "_\\1", rownames(metasample_df_combined)))
  mean_results_list <- list()
  
  # Loop over each unique sample prefix.
  for (sample_prefix in sample_prefixes) {
    # Identify rows corresponding to replicates of the current sample prefix.
    matching_rows <- grep(
      paste0("^", gsub("_auto|_B", "", sample_prefix), "_rep[12]_", 
             sub(".*_", "", sample_prefix)),
      rownames(metasample_df_combined)
    )
    
    # Debugging output: Print sample prefix and its corresponding row names.
    cat("Sample Prefix:", sample_prefix, "\n")
    cat("Matching Rows:\n")
    print(rownames(metasample_df_combined)[matching_rows])
    
    # Merge replicates: If multiple replicates exist, compute their mean.
    if (length(matching_rows) > 1) {
      mean_results_list[[sample_prefix]] <- colMeans(metasample_df_combined[matching_rows, , drop = FALSE])
    } else if (length(matching_rows) == 1) {
      # Use the single replicate directly.
      mean_results_list[[sample_prefix]] <- metasample_df_combined[matching_rows, , drop = FALSE]
    }
  }
  
  # Combine the merged rows into a single data frame.
  if (length(mean_results_list) > 0) {
    merged_df <- do.call(rbind, mean_results_list)
  } else {
    merged_df <- NULL
  }
  return(merged_df)
}

# ------------------------------------------------------------------------------
# Function: plot_metasample_heatmaps
# Description: Creates heatmaps for autologous, monoculture, and combined datasets.
#              For the combined data, both row-split and column-split heatmaps are produced.
#
# Arguments:
#   auto     : A data frame representing the autologous culture metasample matrix.
#   mono     : A data frame representing the monoculture metasample matrix.
#   combined : A data frame representing the combined metasample matrix.
#
# Returns:
#   A list containing four heatmap objects:
#     - heatmap_auto: Heatmap for autologous culture.
#     - heatmap_mono: Heatmap for monoculture.
#     - heatmap_combined_rows: Heatmap for combined data with rows split by culture type.
#     - heatmap_combined_cols: Heatmap for combined data with columns split by culture type.
# ------------------------------------------------------------------------------
plot_metasample_heatmaps <- function(auto, mono, combined) {
  # Convert data frames to matrices.
  auto     <- as.matrix(auto)
  mono     <- as.matrix(mono)
  combined <- as.matrix(combined)
  
  # Define a color palette for the heatmaps.
  heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick"))(100)
  
  # Base arguments for the heatmap generation.
  base_args <- list(
    col = heatmap_colors,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    show_row_names = TRUE,
    show_column_names = TRUE
  )
  
  # Generate heatmap for Autologous culture.
  heatmap_auto <- do.call(Heatmap, c(list(auto, name = "Expression"), base_args))
  
  # Generate heatmap for Monoculture.
  heatmap_mono <- do.call(Heatmap, c(list(mono, name = "Expression"), base_args))
  
  # Generate heatmap for combined data with rows split.
  n_rows <- nrow(combined)
  row_split <- c(rep("Autologous culture", floor(n_rows / 2)),
                 rep("Monoculture", ceiling(n_rows / 2)))
  heatmap_combined_rows <- do.call(Heatmap, c(
    list(combined, name = "Contribution"),
    base_args,
    list(
      cluster_rows = FALSE,
      row_names_gp = gpar(fontsize = 6),
      row_split = row_split,
      row_gap = unit(2, "mm")
    )
  ))
  
  # Generate heatmap for combined data with samples (columns) split.
  combined_t <- t(combined)
  n_cols <- ncol(combined_t)
  column_split <- c(rep("Autologous culture", floor(n_cols / 2)),
                    rep("Monoculture", ceiling(n_cols / 2)))
  heatmap_combined_cols <- do.call(Heatmap, c(
    list(combined_t, name = "Contribution"),
    base_args,
    list(
      cluster_columns = FALSE,
      column_names_gp = gpar(fontsize = 12),
      column_split = column_split,
      column_gap = unit(1, "mm")
    )
  ))
  
  # Return a list of heatmap objects.
  list(
    heatmap_auto = heatmap_auto,
    heatmap_mono = heatmap_mono,
    heatmap_combined_rows = heatmap_combined_rows,
    heatmap_combined_cols = heatmap_combined_cols
  )
}

# ------------------------------------------------------------------------------
# Function: create_filtered_dfs
# Description: Splits a combined data frame into separate data frames for each 
#              culture and patient based on specific patterns in the row names.
#
# Arguments:
#   df: A data frame containing the combined metasample data.
#
# Returns:
#   A list of data frames filtered by culture (autologous/monoculture) and patient (p1, p2, p3).
# ------------------------------------------------------------------------------
create_filtered_dfs <- function(df) {
  df <- as.data.frame(df)
  sn <- rownames(df)
  list(
    auto_p1 = df[grep("_p1_.*_auto", sn), , drop = FALSE],
    auto_p2 = df[grep("_p2_.*_auto", sn), , drop = FALSE],
    auto_p3 = df[grep("_p3_.*_auto", sn), , drop = FALSE],
    mono_p1 = df[grep("_p1_.*_B", sn), , drop = FALSE],
    mono_p2 = df[grep("_p2_.*_B", sn), , drop = FALSE],
    mono_p3 = df[grep("_p3_.*_B", sn), , drop = FALSE]
  )
}

# ------------------------------------------------------------------------------
# Function: process_data
# Description: Processes a data frame by:
#              - Updating row names to extract replicate information.
#              - Extracting time points from row names.
#              - Converting wide-format data to long-format for ggplot2 plotting.
#
# Arguments:
#   data   : A data frame to process (usually one of the filtered data frames).
#   culture: A character string indicating the culture type (e.g., "Autologous" or "Monoculture").
#   patient: A character string indicating the patient identifier.
#
# Returns:
#   A data frame in long format containing the time point, replicate, IC (independent component),
#   value, and metadata (culture and patient) for plotting.
# ------------------------------------------------------------------------------
process_data <- function(data, culture, patient) {
  data <- as.data.frame(data)
  rn <- rownames(data)
  
  # Update row names to retain only the replicate information.
  rownames(data) <- gsub("(rep[12]).*", "\\1", rn)
  
  # Extract the time point from the original row names.
  data$time_point <- gsub("_(p[0-9]+)_.*", "", rn)
  
  # Extract replicate numbers; default to 1 if not found.
  rep_match <- regmatches(rn, regexpr("(?<=_rep)[0-9]+", rn, perl = TRUE))
  data$replicate <- ifelse(nchar(rep_match) > 0, as.numeric(rep_match), 1)
  
  # Convert the data frame from wide to long format.
  pivot_longer(data, cols = starts_with("IC"), names_to = "IC", values_to = "value") %>%
    mutate(
      time_point = factor(time_point, levels = unique(time_point)),
      culture = culture,
      patient = patient
    )
}

# ------------------------------------------------------------------------------
# Function: process_and_plot
# Description: Processes a single group (combination of culture and patient) and generates a
#              line plot showing IC values over time for each replicate.
#
# Arguments:
#   data        : A data frame corresponding to a specific culture and patient.
#   title_prefix: A character string in the form "Autologous - Patient 1" or "Mono - Patient 1",
#                 used to set plot titles and metadata.
#
# Returns:
#   A ggplot2 object representing the line plot of IC values over time for the given group.
# ------------------------------------------------------------------------------
process_and_plot <- function(data, title_prefix) {
  # Parse culture and patient information from the title prefix.
  parts <- unlist(strsplit(title_prefix, " - "))
  culture_label <- ifelse(grepl("Auto", parts[1]), "Autologous", "Monoculture")
  patient_label <- parts[2]
  
  # Process the data into long format including additional metadata.
  long_data <- process_data(data, culture = culture_label, patient = patient_label)
  
  # Generate and return the line plot.
  ggplot(long_data, aes(x = time_point, y = value, group = replicate, color = factor(replicate))) +
    geom_line() +
    facet_wrap(~IC, scales = "fixed", ncol = 5) +
    labs(
      title = paste0(title_prefix, " - IC values over time points"),
      x = "Time point",
      y = "IC value",
      color = "Replicate"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ------------------------------------------------------------------------------
# Function: plot_all_patients_combined
# Description: Merges filtered data frames for all patients and both cultures,
#              then generates a combined line plot comparing IC values over time.
#
# Arguments:
#   filtered_dfs: A list of data frames produced by the create_filtered_dfs() function.
#
# Returns:
#   A ggplot2 object representing the combined line plot for all patients and cultures.
# ------------------------------------------------------------------------------
plot_all_patients_combined <- function(filtered_dfs) {
  combined <- bind_rows(lapply(names(filtered_dfs), function(name) {
    parts <- strsplit(name, "_")[[1]]
    culture <- ifelse(parts[1] == "auto", "Autologous", "Monoculture")
    patient <- paste("Patient", toupper(parts[2]))
    process_data(filtered_dfs[[name]], culture, patient)
  })) %>% arrange(time_point)
  
  ggplot(combined, aes(
    x = time_point, 
    y = value,
    group = interaction(replicate, culture, patient),
    color = culture, 
    linetype = patient
  )) +
    geom_line() +
    facet_wrap(~IC, scales = "fixed", ncol = 5) +
    labs(
      title = "IC values over time points - Autologous vs Monoculture",
      x = "Time point",
      y = "IC value",
      color = "Culture",
      linetype = "Patient"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    scale_linetype_manual(values = c("Patient P1" = "solid",
                                     "Patient P2" = "dashed",
                                     "Patient P3" = "dotted")) +
    scale_color_manual(values = c("Autologous" = "#44AA99",
                                  "Monoculture" = "#E69F00"))
}
