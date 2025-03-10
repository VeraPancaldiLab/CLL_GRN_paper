## ---------------------------------------------------------------------------------------------------------------------------------------
library(OmnipathR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(viper)
library(decoupleR)


## ---------------------------------------------------------------------------------------------------------------------------------------
# Function to get and store database version
get_database_version <- function(database) {
  if (database == "Dorothea" || database == "COLLECTRI") {
    return(packageVersion("decoupleR"))  # decoupleR manages both databases
  } else {
    return("Unknown Database")
  }
}

# Function to save regulons locally (first-time download)
save_database_locally <- function(database, save_path) {
  if (database == "Dorothea") {
    regulons <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))
  } else if (database == "COLLECTRI") {
    regulons <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)
  } else {
    stop("Unsupported database")
  }
  
  # Ensure directory exists
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Save as RDS
  saveRDS(regulons, file = file.path(save_path, paste0(database, "_regulons.rds")))
  return(regulons)
}

# Function to load local database copy
load_database_locally <- function(database, save_path) {
  local_file <- file.path(save_path, paste0(database, "_regulons.rds"))
  if (file.exists(local_file)) {
    regulons <- readRDS(local_file)
  } else {
    stop("Local database file not found. Please run `save_database_locally()` first.")
  }
  return(regulons)
}

# Convert regulons to VIPER format
regulons_to_viper_format <- function(df) {
  if (!"source" %in% colnames(df)) {
    stop("Error: 'source' column is missing from the regulons data.")
  }
  regulon_list <- split(df, df$source)
  viper_regulons <- lapply(regulon_list, function(regulon) {
    tfmode <- stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  return(viper_regulons)
}

# Function to compute TF scores with VIPER
compute_TFs_scores <- function(RNA.counts.normalized, output_dir, database, save_path) {
  
  # Load the local version of the database
  regulons <- load_database_locally(database, save_path)
  
  # Convert regulons to VIPER format
  regu <- regulons_to_viper_format(regulons)
  
  # Convert RNA counts to a data frame (if not already in that format)
  RNA.counts.normalized <- as.data.frame(RNA.counts.normalized)
  
  # Run VIPER
  vpres <- viper(RNA.counts.normalized, regu, verbose = FALSE, minsize = 4)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the VIPER results
  file_suffix <- ifelse(database == "Dorothea", "TF_score_dorothea", "TF_score_collectri")
  output_path <- file.path(output_dir, paste0(file_suffix, ".csv"))
  write.csv(vpres, output_path)
  
  return(vpres)
}
