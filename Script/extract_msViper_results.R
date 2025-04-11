# Auto resutls
library(tidyverse)

# Function to read all RDS files from a folder
read_all_rds <- function(folder_path) {
  # Get list of .rds files in the folder
  rds_files <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)
  
  # Read each file and store in a named list
  rds_data <- lapply(rds_files, readRDS)
  names(rds_data) <- basename(rds_files) %>% str_remove("\\.rds$")
  
  return(rds_data)
}

# Specify the folder path
folder_path <- "Script/msViper/Results/Monoculture/"

# Call the function to read all RDS files
all_rds_data <- read_all_rds(folder_path)

filtered_results <- list()
list_TFs_all_patients <- c()

# Loop through each object in the all_rds_data list
for (file_name in names(all_rds_data)) {
  object <- all_rds_data[[file_name]]
  
    # Identify elements where p.value < 0.05
    significant_indices <- which(object$p.value < 0.05)
    
    significant_TFs <- names(object$p.value)[significant_indices]
    
    # Determine the source from the file name (e.g., 1vs4, 4vs8, etc.)
    source <- str_extract(file_name, "1vs4|4vs8|8vs11|11vs14")
    patient <- str_extract(file_name, "Patient1|Patient2|Patient3")
    
    # Extract significant elements
    filtered_data <- list(
      nes = object$nes[significant_indices],
      size = object$size[significant_indices],
      p.value = object$p.value[significant_indices]
    )
    
    # Add filtered data to the results list
    filtered_results[[file_name]] <- filtered_data
   
     # Create a dataframe for this file's results
    file_TFs_df <- data.frame(TF = significant_TFs, Source = source, Patient = patient, stringsAsFactors = FALSE)
    
    # Append to the main dataframe
    list_TFs_all_patients <- rbind(list_TFs_all_patients, file_TFs_df)

}
auto_TFs <- list_TFs_all_patients
write.csv(auto_TFs, "Results/Temp/msViper_auto_TFs_comparison_info.csv")
auto_TFs_msViper <- unique(auto_TFs$TF)
write.csv(auto_TFs_msViper, "Results/Temp/msViper_auto_TFs.csv")

mono_TFs <- list_TFs_all_patients
write.csv(mono_TFs, "Results/Temp/msViper_mono_TFs_comparison_info.csv")
mono_TFs_msViper <- unique(mono_TFs$TF)
write.csv(mono_TFs_msViper, "Results/Temp/msViper_mono_TFs.csv")
# Run again the code above for the two conditions and take the merge between the two lists 

all_TFs_msViper <- unique(c(auto_TFs_msViper, mono_TFs_msViper))

Patient1_TFs <- rbind((auto_TFs %>% filter(auto_TFs$Patient == "Patient1")), (mono_TFs %>% filter(mono_TFs$Patient == "Patient1")))
Patient1_TFs <- unique(Patient1_TFs$TF)
Patient2_TFs <- rbind((auto_TFs %>% filter(auto_TFs$Patient == "Patient2")), (mono_TFs %>% filter(mono_TFs$Patient == "Patient2")))
Patient2_TFs <- unique(Patient2_TFs$TF)
Patient3_TFs <- rbind((auto_TFs %>% filter(auto_TFs$Patient == "Patient3")), (mono_TFs %>% filter(mono_TFs$Patient == "Patient3")))
Patient3_TFs <- unique(Patient3_TFs$TF)
# msViper_TFs <- unique(c(auto_TFs_msViper, mono_TFs_msViper))

# Load the TF activity matrices for auto and mono
TFact_auto <- read.csv("Results/Viper/Autologous/TF_score_collectri.csv", row.names = 1) 
TFact_mono <- read.csv("Results/Viper/B_cell/TF_score_collectri.csv", row.names = 1)

# Filter these matrices to select only the most DA TFs and then do the merge
TFact_auto_DA_TFs <- TFact_auto[most_DA_TFs$`unique(c(auto_TFs$TF, mono_TFs$TF))`, , drop = FALSE]
TFact_mono_DA_TFs <- TFact_mono[most_DA_TFs$`unique(c(auto_TFs$TF, mono_TFs$TF))`, , drop = FALSE]

TFact_DA_TFs <- cbind(TFact_auto_DA_TFs, TFact_mono_DA_TFs)
# Extract unique time points, patients, and conditions
time_points <- c("D1", "D4", "D8", "D11", "D14")
patients <- c("p1", "p2", "p3")
conditions <- c("auto", "B")

# Initialize an empty dataframe for the means
mean_df <- data.frame(row_id = seq_len(nrow(df)))
rownames(mean_df) <- rownames(df)

# Loop through time points, patients, and conditions to compute means
for (time in time_points) {
  for (patient in patients) {
    for (condition in conditions) {
      # Build column names for replicates
      rep1_col <- paste0(time, "_", patient, "_rep1_", condition)
      rep2_col <- paste0(time, "_", patient, "_rep2_", condition)
      
      # Check if both columns exist in the dataframe
      if (rep1_col %in% colnames(df) && rep2_col %in% colnames(df)) {
        # Compute mean of replicates
        mean_col <- rowMeans(df[, c(rep1_col, rep2_col)])
        
        # Add the mean to the new dataframe
        mean_df[[paste0(time, "_", patient, "_mean_", condition)]] <- mean_col
      }
    }
  }
}

# Remove the row_id column (if not needed)
mean_df$row_id <- NULL

correlation_matrix <- cor(t(mean_df)) %>% melt() %>% filter(abs(value) >= 0.5)

#######################################################################

# Compare with the TFs in the GRN
TFs_previous_GRN <- read.csv("Results/Temp/TFs_auto_mono_separately.csv")

common_Tfs <- intersect(msViper_TFs, TFs_previous_GRN$x)
