# ===============================================================================
# BATCH PROCESSING SCRIPT FOR VPD DROUGHT ANALYSIS - MULTIPLE FLUXNET SITES
# ===============================================================================
# This script processes multiple FLUXNET sites for VPD drought analysis
# It reads the folder list and processes each site sequentially
# ===============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries
library(tidyverse)
library(lubridate)

# ===============================================================================
# BATCH CONFIGURATION
# ===============================================================================

BATCH_CONFIG <- list(
  # Directory containing all site-specific ERA5/ERA5-Land input folders.
  # Please replace this with your local data directory.
  #
  # Example folder structure:
  # ~/data/atmospheric_drought_inputs/
  # ├── DK_Sor_1996_2020/
  # │   └── DK-Sor_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv
  # ├── FR_Pue_2000_2020/
  # │   └── FR-Pue_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv
  # └── ...
  base_dir = "~/data/atmospheric_drought_inputs",
  
  # CSV file containing the list of site folders to process.
  # The CSV should contain a column named "folder".
  #
  # Example:
  # folder
  # DK_Sor_1996_2020
  # FR_Pue_2000_2020
  folder_list_file = "folder_list.csv",
  
  # Directory where atmospheric drought output files will be saved.
  # Please replace this with your preferred output directory.
  output_base_dir = "~/data/atmospheric_drought_outputs",
  
  # Log file for batch processing.
  log_file = "atmospheric_drought_batch_processing_log.txt",
  
  # Source file containing the atmospheric drought analysis functions.
  # Use the single-site example script, but the batch script will skip
  # the example execution line at the bottom.
  source_file = "AtmosphericDrought_single_site_example.R"
)

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

# Function to extract site code from folder name
extract_site_code <- function(folder_name) {
  # Extract first two characters (country code)
  country_code <- substr(folder_name, 1, 2)
  
  # Find the position of the underscore
  underscore_pos <- regexpr("_", folder_name)[1]
  
  # Extract site code (3 characters after underscore)
  if (underscore_pos > 0) {
    site_code <- substr(folder_name, underscore_pos + 1, underscore_pos + 3)
  } else {
    stop("No underscore found in folder name")
  }
  
  # Combine with hyphen
  full_code <- paste0(country_code, "-", site_code)
  
  return(list(
    country_code = country_code,
    site_code = site_code,
    full_code = full_code
  ))
}

# Function to process a single site for VPD drought
process_single_site_vpd <- function(folder_name, batch_config, log_con) {
  
  # Extract site codes
  site_codes <- extract_site_code(folder_name)
  
  cat("\n========================================\n")
  cat("Processing VPD Drought Analysis for:", folder_name, "\n")
  cat("Site code:", site_codes$full_code, "\n")
  cat("========================================\n")
  
  writeLines(paste("\n=== Processing VPD Drought for", folder_name, "==="), log_con)
  
  # Create site-specific configuration
  MASTER_CONFIG <- list(
    # Working directory - the specific site folder
    working_dir = file.path(batch_config$base_dir, folder_name),
    
    # Input data file - constructed based on site code
    input_file = paste0(site_codes$full_code, "_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv"),
    
    # Moving window size (days on each side of target day)
    window_size = 15,
    
    # Quantile for threshold calculation (95th percentile for high VPD)
    threshold_quantile = 0.95,
    
    # Output file prefix - save to results directory
    output_prefix = file.path(batch_config$output_base_dir, folder_name, paste0(site_codes$full_code, "_VPD_Drought")),
    
    # Vine copula family set
    copula_families = c(1:10),
    
    # Selection criterion for vine copula
    selection_criterion = "AIC"
  )
  
  # Period configurations (same for all sites) - UPDATED
  PERIOD_CONFIGS <- list(
    # Period 1: Analysis for 1991-2000 based on 1961-1990 reference
    period_1 = list(
      reference_start = 1961,
      reference_end = 1990,
      target_start = 1991,
      target_end = 2000
    ),
    # Period 2: Analysis for 2001-2010 based on 1971-2000 reference
    period_2 = list(
      reference_start = 1971,
      reference_end = 2000,
      target_start = 2001,
      target_end = 2010
    ),
    # Period 3: Analysis for 2011-2020 based on 1981-2010 reference
    period_3 = list(
      reference_start = 1981,
      reference_end = 2010,
      target_start = 2011,
      target_end = 2020
    ),
    # Period 4: Analysis for 2021-2025 based on 1990-2020 reference
    period_4 = list(
      reference_start = 1990,
      reference_end = 2020,
      target_start = 2021,
      target_end = 2024
    )
  )
  
  # Create output directory for this site
  site_output_dir <- file.path(batch_config$output_base_dir, folder_name)
  if (!dir.exists(site_output_dir)) {
    dir.create(site_output_dir, recursive = TRUE)
  }
  
  # Try to process the site
  tryCatch({
    # Check if input file exists
    input_path <- file.path(MASTER_CONFIG$working_dir, MASTER_CONFIG$input_file)
    if (!file.exists(input_path)) {
      stop(paste("Input file not found:", MASTER_CONFIG$input_file))
    }
    
    # Check if ERA5 data contains VPD_kPa column
    test_data <- read.csv(input_path, nrows = 10)
    if (!"VPD_kPa" %in% names(test_data)) {
      stop("VPD_kPa column not found in ERA5 data file")
    }
    
    # Run the VPD drought analysis
    results <- run_combined_vpd_drought_analysis(MASTER_CONFIG, PERIOD_CONFIGS)
    
    writeLines(paste("  SUCCESS: VPD Drought analysis completed for", folder_name), log_con)
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    writeLines(paste("  ERROR:", e$message), log_con)
    return(FALSE)
  })
}

# ===============================================================================
# MAIN EXECUTION
# ===============================================================================

# Create output directory if it doesn't exist
if (!dir.exists(BATCH_CONFIG$output_base_dir)) {
  dir.create(BATCH_CONFIG$output_base_dir, recursive = TRUE)
}

# Set working directory to base directory
setwd(BATCH_CONFIG$base_dir)

# Source the main VPD drought analysis functions
# First, we need to prevent the original script from running its example
# We'll do this by defining a flag
BATCH_MODE <- TRUE

# Now source the file, but skip the execution part
source_lines <- readLines(BATCH_CONFIG$source_file)

# Find the line that starts the execution (look for the run_combined_vpd_drought_analysis call)
exec_line <- which(grepl("^drought_results <- run_combined_vpd_drought_analysis", source_lines))
if (length(exec_line) > 0) {
  # Source only up to before the execution
  source_text <- paste(source_lines[1:(exec_line-1)], collapse = "\n")
  eval(parse(text = source_text))
} else {
  # If we can't find the execution line, source the whole file
  source(BATCH_CONFIG$source_file)
}

# Read folder list
folder_list <- read.csv(BATCH_CONFIG$folder_list_file)
cat("Found", nrow(folder_list), "sites to process for VPD Drought\n\n")

# Open log file
log_con <- file(BATCH_CONFIG$log_file, open = "wt")
writeLines(paste("VPD Drought batch processing started at:", Sys.time()), log_con)
writeLines(paste("Total sites:", nrow(folder_list)), log_con)

# Initialize counters
successful <- 0
failed <- 0
skipped <- 0

# Process each site
for (i in 1:nrow(folder_list)) {
  folder_name <- folder_list$folder[i]
  
  cat("\n###############################################\n")
  cat("VPD Drought Analysis - Site", i, "of", nrow(folder_list), "\n")
  
  # Check if we've already processed this site
  site_codes <- extract_site_code(folder_name)
  expected_output <- file.path(BATCH_CONFIG$output_base_dir, folder_name, 
                               paste0(site_codes$full_code, "_VPD_Drought_1991_2024.csv"))
  
  if (file.exists(expected_output)) {
    cat("Site already processed for VPD drought, skipping:", folder_name, "\n")
    writeLines(paste("\nSkipping", folder_name, "(already processed)"), log_con)
    skipped <- skipped + 1
    next
  }
  
  # Process the site
  success <- process_single_site_vpd(folder_name, BATCH_CONFIG, log_con)
  
  if (success) {
    successful <- successful + 1
  } else {
    failed <- failed + 1
  }
  
  # Flush log
  flush(log_con)
  
  # Optional: Add a small delay between sites
  Sys.sleep(1)
}

# Final summary
cat("\n\n=======================================\n")
cat("VPD DROUGHT BATCH PROCESSING COMPLETE\n")
cat("=======================================\n")
cat("Total sites:", nrow(folder_list), "\n")
cat("Successful:", successful, "\n")
cat("Failed:", failed, "\n")
cat("Skipped:", skipped, "\n")

writeLines("\n=======================================", log_con)
writeLines("VPD DROUGHT BATCH PROCESSING SUMMARY", log_con)
writeLines(paste("Completed at:", Sys.time()), log_con)
writeLines(paste("Successful:", successful), log_con)
writeLines(paste("Failed:", failed), log_con)
writeLines(paste("Skipped:", skipped), log_con)

# Close log file
close(log_con)

cat("\nCheck", BATCH_CONFIG$log_file, "for details.\n")

# Create a summary CSV of the batch processing
summary_df <- data.frame(
  folder_name = folder_list$folder,
  processed = NA,
  status = NA,
  has_vpd_data = NA,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(summary_df)) {
  folder_name <- summary_df$folder_name[i]
  site_codes <- extract_site_code(folder_name)
  expected_output <- file.path(BATCH_CONFIG$output_base_dir, folder_name, 
                               paste0(site_codes$full_code, "_VPD_Drought_1991_2024.csv"))
  
  if (file.exists(expected_output)) {
    summary_df$processed[i] <- TRUE
    summary_df$status[i] <- "Success"
    summary_df$has_vpd_data[i] <- TRUE
  } else {
    # Check if ERA5 file exists and has VPD data
    era5_file <- file.path(BATCH_CONFIG$base_dir, folder_name, 
                           paste0(site_codes$full_code, "_ERA5_Land_stitched_1970-12-14_to_2022-12-31.csv"))
    if (file.exists(era5_file)) {
      # Quick check for VPD column
      tryCatch({
        test_data <- read.csv(era5_file, nrows = 5)
        summary_df$has_vpd_data[i] <- "VPD_kPa" %in% names(test_data)
      }, error = function(e) {
        summary_df$has_vpd_data[i] <- FALSE
      })
    } else {
      summary_df$has_vpd_data[i] <- FALSE
    }
    
    summary_df$processed[i] <- FALSE
    summary_df$status[i] <- ifelse(summary_df$has_vpd_data[i], "Failed", "No VPD Data")
  }
}

# Save summary
summary_file <- file.path(BATCH_CONFIG$output_base_dir, "vpd_drought_batch_processing_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat("\nVPD Drought processing summary saved to:", summary_file, "\n")

# Display sites that couldn't be processed due to missing VPD data
no_vpd_sites <- summary_df[!summary_df$has_vpd_data,]
if (nrow(no_vpd_sites) > 0) {
  cat("\n=======================================\n")
  cat("Sites without VPD data:\n")
  print(no_vpd_sites[, c("folder_name", "status")])
  cat("\nThese sites need VPD_kPa column in their ERA5 data files.\n")
}
