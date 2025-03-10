## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: load_msViper_files
# Description:
#   Reads a set of msVIPER result files (in RDS format), filters the data based on a p-value threshold,
#   extracts transcription factor names (from row names) and their normalized enrichment scores (NES).
#
# Arguments:
#   files             - A character vector of file paths to the msVIPER RDS files.
#   p_value_threshold - Numeric; threshold for filtering rows by the 'p.value' column (default is 0.05).
#
# Returns:
#   A list with three elements:
#     data - A list of data frames, one per file, with rows meeting the p-value criterion.
#     tfs  - A list of transcription factor names (from rownames) for each file.
#     nes  - A list of NES values (with names set to the TF names) for each file.

load_msViper_files <- function(files, p_value_threshold = 0.05) {
  # Initialize empty lists to store the data, transcription factor names, and NES values.
  d_list <- list()
  tfs_list <- list()
  nes_list <- list()
  
  # Loop over each file provided in the files vector.
  for (f in files) {
    # Read the RDS file and convert it to a data frame.
    # The code assumes the second element is not needed, so it's removed with [-2].
    ddf <- as.data.frame(readRDS(f)[-2])
    
    # Filter the data frame to keep only rows where the p.value is below the threshold.
    d_list[[f]] <- ddf[which(ddf[,'p.value'] < p_value_threshold), ]
    
    # Extract transcription factor names from the row names of the filtered data frame.
    tfs_list[[f]] <- rownames(d_list[[f]])
    
    # Extract the normalized enrichment scores (NES) from the 'nes' column.
    nes_list[[f]] <- d_list[[f]]$nes
    
    # Set the names of the NES vector to be the transcription factor names.
    names(nes_list[[f]]) <- tfs_list[[f]]
  }
  
  # Return a list containing the filtered data, transcription factor names, and NES values.
  return(list(data = d_list, tfs = tfs_list, nes = nes_list))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: build_tf_matrix
# Description:
#   Constructs a transcription factor (TF) matrix from msVIPER results for both autologous and monoculture conditions.
#   The function aligns the TF names and their corresponding NES values, and calculates a correlation matrix across conditions.
#
# Arguments:
#   tfs_auto   - List of transcription factor names for the autologous condition.
#   tfs_mono   - List of transcription factor names for the monoculture condition.
#   nes_auto   - List of NES values for the autologous condition.
#   nes_mono   - List of NES values for the monoculture condition.
#   files_auto - Character vector of file paths for the autologous condition.
#   files_mono - Character vector of file paths for the monoculture condition.
#
# Returns:
#   A list containing:
#     mat           - Matrix with rows as unique TFs and columns as each sample/file.
#     allcondcor    - Correlation matrix computed across conditions (columns).
#     tfsuniqueauto - Unique TF names identified from the autologous data.
#     tfsuniquemono - Unique TF names identified from the monoculture data.
#     tfsunique     - Combined unique TF names from both conditions.
#     nes           - Combined NES values.
#     tfs           - Combined TF names.
build_tf_matrix <- function(tfs_auto, tfs_mono, nes_auto, nes_mono, files_auto, files_mono) {
  # Combine the unique TF names for the autologous condition
  tfsuniqueauto <- unique(unlist(tfs_auto[files_auto]))
  # Combine the unique TF names for the monoculture condition
  tfsuniquemono <- unique(unlist(tfs_mono[files_mono]))
  
  # Create a combined unique list of TF names from both conditions
  tfsunique <- unique(unlist(list(auto = tfsuniqueauto, mono = tfsuniquemono)))
  
  # Combine the NES and TF lists from both conditions into single lists
  nes <- c(nes_auto, nes_mono)
  tfs <- c(tfs_auto, tfs_mono)
  
  # Combine all file names (order here is important if you need a particular order later)
  files <- c(files_mono, files_auto)
  
  # Create an empty matrix with rows as unique TFs and columns as each file
  mat <- matrix(0, nrow = length(tfsunique), ncol = length(files))
  colnames(mat) <- files
  rownames(mat) <- tfsunique
  
  # Loop through each file and assign the corresponding NES values into the matrix
  # using the TF names as row indices
  for (f in files) {
    # 'tfs[[f]]' should contain the row names (TFs) for this file
    # and 'nes[[f]]' holds the corresponding NES values.
    mat[tfs[[f]], f] <- nes[[f]]
  }
  
  # Rename the columns so that files with "Monoculture" are prefixed with "Mono_" and 
  # files with autologous data with "Auto_". The regex extracts the patient and time point.
  colnames(mat) <- paste0(
    ifelse(grepl("Monoculture", colnames(mat)), "Mono_", "Auto_"),
    sub(".*msViper_result_(Patient\\d+)_TimePoint_(\\d+vs\\d+).*", "\\1_\\2", colnames(mat))
  )
  
  # Compute a correlation matrix of the transposed matrix (correlation among conditions)
  allcondcor <- cor(t(mat))
  
  # Optionally reorder columns using natural (mixed) order if the gtools package is available
  if (requireNamespace("gtools", quietly = TRUE)) {
    mat <- mat[, gtools::mixedorder(colnames(mat))]
  }
  
  # Return both the matrix and its correlation matrix as a list
  return(list(mat = mat, allcondcor = allcondcor, tfsuniqueauto = tfsuniqueauto, tfsuniquemono = tfsuniquemono, tfsunique = tfsunique, nes = nes, tfs = tfs))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function: calculate_correlation_and_export
# Description:
#   Calculates correlation matrices based on the TF activity matrix, generates a network graph
#   from TF correlations, and exports both the edge list of the graph and the TF annotation matrix to files.
#
# Arguments:
#   mat    - A matrix where rows are TFs and columns are samples/files with NES values.
#   outdir - (Optional) Output directory for saving the exported files (default is "Results/msVIPER").
#
# Returns:
#   A list containing:
#     corcond - Correlation matrix between conditions (columns of the matrix).
#     cortf   - Correlation matrix between transcription factors (rows of the matrix) computed by transposing the matrix.
#     graph   - An igraph object representing the undirected weighted network graph of TF correlations.
#     mat     - The original TF annotation matrix used for the analysis.
calculate_correlation_and_export <- function(mat, outdir = "Results/msVIPER") {
  # Calculate the correlation between conditions (columns)
  corcond <- cor(mat)
  
  # Calculate the correlation between TFs (rows) by transposing the matrix
  cortf <- cor(t(mat))
  
  # Create an undirected weighted graph from the TF correlation matrix using igraph
  Gcortf <- igraph::graph_from_adjacency_matrix(cortf, weight = TRUE, mode = "undirected")
  
  # Write the edge list (with weights)
  write.table(cbind(as_edgelist(Gcortf), igraph::E(Gcortf)$weight),
              file = file.path(outdir, "Gcortfnet.txt"),
              quote = FALSE, row.names = FALSE)
  
  # Write the TF annotation matrix
  write.table(mat, file = file.path(outdir, "TFannot.txt"), 
              quote = FALSE, sep = "\t")
  
  # Write the correlation among TF activities
  write.table(cortf, file = file.path(outdir, "TFcorr.txt"), 
              quote = FALSE, sep = "\t")
  
  # Optionally, return the correlation matrices and graph for further inspection
  return(list(corcond = corcond, cortf = cortf, graph = Gcortf, mat = mat))
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function to process culture data (either autologous or monoculture)
# Arguments:
#   files: a character vector of file paths (already sorted if needed)
#   tfs_list: a named list where each element (named by file) is a vector of TF names from that file
#   nes_list: a named list where each element (named by file) is a vector of NES values (names must match TF names)
#   tfs_unique: a vector of unique TF names to use as row names in the matrix
# Returns:
#   A list containing:
#     full: the full TF matrix (rows = TFs, columns = files)
#     patient1, patient2, patient3: submatrices for each patient
process_culture_data <- function(files, tfs_list, nes_list, tfs_unique) {
  
  # Order files in natural (mixed) order.
  files <- mixedsort(files)
  
  # Initialize a matrix with rows as unique TFs and columns corresponding to the files.
  mat <- matrix(0, nrow = length(tfs_unique), ncol = length(files))
  colnames(mat) <- files
  rownames(mat) <- tfs_unique
  
  # Fill the matrix with NES values. For each file, assign NES values to rows given by its TF names.
  for (f in files) {
    mat[tfs_list[[f]], f] <- nes_list[[f]]
  }
  
  # Rename the columns using a regex that extracts patient and time information.
  # It prefixes with "Mono_" if "Monoculture" is in the file name; otherwise "Auto_".
  colnames(mat) <- paste0(
    ifelse(grepl("Monoculture", colnames(mat)), "Mono_", "Auto_"),
    sub(".*msViper_result_(Patient\\d+)_TimePoint_(\\d+vs\\d+).*", "\\1_\\2", colnames(mat))
  )
  
  # Reorder the columns using natural (mixed) order.
  mat <- mat[, mixedorder(colnames(mat))]
  
  # Extract submatrices for each patient by selecting columns with corresponding patient strings.
  mat_p1 <- mat[, grep("Patient1", colnames(mat))]
  mat_p1 <- mat_p1[which(rowSums(mat_p1) > 0), ]
  
  mat_p2 <- mat[, grep("Patient2", colnames(mat))]
  mat_p2 <- mat_p2[which(rowSums(mat_p2) > 0), ]
  
  mat_p3 <- mat[, grep("Patient3", colnames(mat))]
  mat_p3 <- mat_p3[which(rowSums(mat_p3) > 0), ]
  
  return(list(full = mat, patient1 = mat_p1, patient2 = mat_p2, patient3 = mat_p3))
}

# Function to plot heatmaps for the processed culture data.
# Arguments:
#   data_list: a list with elements 'full', 'patient1', 'patient2', and 'patient3' (matrices)
#   culture_label: a string to label the culture type in the title (e.g., "Autologous" or "Monoculture")
plot_culture_data <- function(data_list, culture_label) {
  
  # Plot full heatmap for all patients
  pheatmap(data_list$full, fontsize_row = 8, cluster_cols = FALSE, 
         main = paste("All patients -", culture_label), 
         scale = "none")
    
  # Plot heatmap for each patient
  pheatmap(data_list$patient1, fontsize_row = 8, cluster_cols = FALSE, 
         main = paste("Patient 1 -", culture_label), 
         scale = "none")
  
  pheatmap(data_list$patient2, fontsize_row = 8, cluster_cols = FALSE, 
         main = paste("Patient 2 -", culture_label), 
         scale = "none")
  
  pheatmap(data_list$patient3, fontsize_row = 8, cluster_cols = FALSE, 
         main = paste("Patient 3 -", culture_label), 
         scale = "none")
}


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function to annotate TF indices for a given patient
# Arguments:
#   full_mat: The full matrix with TFs as rows and conditions as columns.
#   patient_label: A string (e.g., "Patient1") to select columns belonging to that patient.
#
# Returns:
#   A list with two components:
#     positive: A list (one element per column) containing indices of TFs with positive values.
#     negative: A list (one element per column) containing indices of TFs with negative values.
annotate_tfs_by_patient <- function(full_mat, patient_label) {
  # Subset the matrix to columns matching the patient label
  patient_mat <- full_mat[, grep(patient_label, colnames(full_mat))]
  
  # Keep only rows that have at least one non-zero value
  patient_mat <- patient_mat[rowSums(patient_mat) > 0, ]
  
  # For each column, find indices of TFs with positive NES values
  pos_annotation <- apply(patient_mat, 2, function(x) which(x > 0))
  
  # For each column, find indices of TFs with negative NES values
  neg_annotation <- apply(patient_mat, 2, function(x) which(x < 0))
  
  return(list(positive = pos_annotation, negative = neg_annotation))
}

