# =============================================================================
# Carbon Flux Anomaly Analysis
#
# Purpose:
#   This script calculates daily carbon flux anomalies for multiple sites.
#   For each site, climatological seasonal cycles are calculated for NEE, GPP,
#   and ecosystem respiration using a month-day based 31-day moving window.
#   Anomalies are calculated after removing both the climatological seasonal
#   cycle and a linear temporal trend.
#
# Method:
#   1. Apply quality-control filtering using NEE_VUT_REF_QC.
#   2. Calculate climatology using a month-day based 31-day moving window.
#   3. Calculate residuals as:
#        residual = observed value - climatological mean
#   4. Fit and remove a linear trend from the residuals.
#   5. Calculate anomalies and standardized z-scores.
#
# Input data structure:
#   base_dir/
#     Folder_List_All.csv
#     SITE_YYYY_YYYY/
#       SITE.csv
#
# Example:
#   base_dir/
#     Folder_List_All.csv
#     AU_How_2001_2014/
#       AU_How.csv
#
# Required columns in each site CSV:
#   TIMESTAMP
#   NEE_VUT_REF
#   NEE_VUT_REF_QC
#   GPP_NT_VUT_REF
#   RECO_NT_VUT_REF
#
# Main output:
#   carbon_anomalies_QC_filtered_SITE.csv
#
# Notes:
#   - The baseline period for each site is extracted from the folder name,
#     e.g. AU_How_2001_2014 uses 2001-2014 as the baseline period.
#   - The site CSV file is assumed to be named using the first six characters
#     of the folder name, e.g. AU_How_2001_2014/AU_How.csv.
#   - Missing values coded as -9999 are treated as NA.
# =============================================================================


rm(list = ls())

# =============================================================================
# Section 0: Libraries
# =============================================================================

library(tidyverse)
library(lubridate)
library(ggplot2)
library(mgcv)
library(rlang)


# =============================================================================
# Section 1: User Configuration
# =============================================================================

# Path to the folder containing Folder_List_All.csv and all site folders.
# Please replace this with the path to your local data directory.
#
# Example data structure:
#   base_dir/
#     Folder_List_All.csv
#     AU_How_2001_2014/
#       AU_How.csv
#     BE_Bra_1997_2014/
#       BE_Bra.csv

base_dir <- "path/to/your/data_directory"

# CSV file listing all site folders.
# The file must contain a column named "folder".
folder_list_file <- "Folder_List_All.csv"

# Path to the folder where output files will be saved.
# Please replace this with your preferred output directory.
output_base_dir <- "path/to/your/output_directory"

# Log file name.
log_file <- "carbon_anomaly_batch_processing_log.txt"

if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir, recursive = TRUE)
}

setwd(base_dir)
folder_list <- read.csv(folder_list_file)

# Open the log file
log_con <- file(log_file, open = "wt")
writeLines(paste("Carbon Anomaly batch processing started at:", Sys.time()), log_con)
writeLines(paste("Total sites:", nrow(folder_list)), log_con)

# Initialize counters
successful <- 0
failed <- 0
skipped <- 0

# Process each site
for (i in 1:nrow(folder_list)) {
  folder_name <- folder_list$folder[i]
  
  # Extract the first six characters as the file prefix
  file_prefix <- substr(folder_name, 1, 6)
  
  # Extract the start and end years
  start_year <- as.integer(substr(folder_name, 8, 11))
  end_year <- as.integer(substr(folder_name, 13, 16))
  
  cat("\nSite", i, "/", nrow(folder_list), ":", folder_name, "\n")
  
  # Check whether the site has already been processed
  expected_output <- file.path(output_base_dir, folder_name, 
                               paste0("carbon_anomalies_QC_filtered_", file_prefix, ".csv"))
  
  if (file.exists(expected_output)) {
    cat("Already processed, skipping\n")
    writeLines(paste("\nSkipping", folder_name, "(already processed)"), log_con)
    skipped <- skipped + 1
    next
  }
  
  writeLines(paste("\n=== Processing", folder_name, "==="), log_con)
  
  site_dir <- file.path(base_dir, folder_name)
  site_output_dir <- file.path(output_base_dir, folder_name)
  
  if (!dir.exists(site_output_dir)) {
    dir.create(site_output_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Set the working directory to the site folder
    setwd(site_dir)
    
    # =============================================================================
    # Start of the month-day workflow
    # =============================================================================
    
    # --- File and Variable Settings ---
    file_path         <- paste0(file_prefix, ".csv")
    target_vars       <- c("NEE_VUT_REF", "GPP_NT_VUT_REF", "RECO_NT_VUT_REF")
    qc_variable       <- "NEE_VUT_REF_QC"
    qc_threshold      <- 0.7
    
    # --- Time Period Settings ---
    baseline_years    <- start_year:end_year
    analysis_years    <- start_year:end_year
    detrend_years     <- start_year:end_year
    export_years      <- start_year:end_year
    
    # --- Method Settings ---
    window_size       <- 31
    extreme_threshold <- 2
    
    # Check whether the window size is an odd number
    if (window_size %% 2 == 0) {
      stop("Error: window_size must be an odd number")
    }
    window_offset <- (window_size - 1) / 2
    
    # Read the data
    if (!file.exists(file_path)) {
      stop("File not found: ", file_path)
    }
    
    data_raw <- read.csv(file_path, stringsAsFactors = FALSE, na.strings = "-9999")
    
    # Check whether all target variables and the QC variable exist
    required_cols <- c(target_vars, qc_variable)
    if (!all(required_cols %in% names(data_raw))) {
      missing_vars <- required_cols[!required_cols %in% names(data_raw)]
      stop(paste("Error: The following columns were not found in the file:", paste(missing_vars, collapse=", ")))
    }
    
    # Keep only the required columns
    data <- data_raw[, c("TIMESTAMP", required_cols)]
    
    # Convert the timestamp to Date format
    data$Date <- as.Date(as.character(data$TIMESTAMP), format = "%Y%m%d")
    data <- data[!is.na(data$Date), ]
    
    # Apply the quality-control filter
    rows_before_qc <- nrow(data)
    data <- data %>% filter(!!sym(qc_variable) >= qc_threshold)
    
    data$Year <- year(data$Date)
    data$Month <- month(data$Date)
    data$Day <- day(data$Date)
    data$DOY <- yday(data$Date)  # Used only for plotting
    
    # Create a month-day string as the key
    data$MonthDay <- sprintf("%02d-%02d", data$Month, data$Day)
    
    # Initialize a data frame to store all results
    all_results_data <- data
    
    # Create a PDF file to save the plots
    pdf(file.path(site_output_dir, paste0(file_prefix, "_carbon_anomaly_plots.pdf")))
    
    # Perform the analysis for each variable
    for (current_var in target_vars) {
      
      # Create a temporary working data frame for the current variable and remove missing values
      working_data <- all_results_data %>%
        select(Date, Year, Month, Day, MonthDay, DOY, all_of(current_var)) %>%
        rename(target_var = !!sym(current_var)) %>%
        filter(!is.na(target_var))
      
      # Calculate climatology using the month-day method
      climatology_data <- working_data %>% filter(Year %in% baseline_years)
      
      # Calculate the climatological mean for each month-day combination
      calculate_monthday_mean <- function(month, day, data_for_calc) {
        target_monthday <- sprintf("%02d-%02d", month, day)
        
        # Get all dates within the moving window
        window_dates <- c()
        
        for (offset in -window_offset:window_offset) {
          window_date <- as.Date(paste(2020, month, day, sep="-")) + offset  # Use a leap year as the reference
          window_month <- month(window_date)
          window_day <- day(window_date)
          window_monthday <- sprintf("%02d-%02d", window_month, window_day)
          
          # Special handling for February 29
          if (window_monthday == "02-29") {
            # Check whether February 29 exists in the data, i.e. whether leap-year data are available
            if (any(data_for_calc$MonthDay == "02-29")) {
              window_dates <- c(window_dates, window_monthday)
            }
            # If February 29 does not exist in the data, skip this day
          } else {
            window_dates <- c(window_dates, window_monthday)
          }
        }
        
        window_subset <- data_for_calc %>% filter(MonthDay %in% window_dates)
        
        if (nrow(window_subset) > 0) {
          return(mean(window_subset$target_var, trim = 0.1, na.rm = TRUE))
        } else {
          return(NA)
        }
      }
      
      # Create all possible month-day combinations
      all_monthdays <- unique(climatology_data$MonthDay)
      monthday_climatology <- data.frame()
      
      for (md in all_monthdays) {
        month_val <- as.integer(substr(md, 1, 2))
        day_val <- as.integer(substr(md, 4, 5))
        mean_val <- calculate_monthday_mean(month_val, day_val, climatology_data)
        monthday_climatology <- rbind(monthday_climatology, 
                                      data.frame(MonthDay = md, 
                                                 climatology_mean = mean_val))
      }
      
      # Calculate residuals
      working_data <- working_data %>% left_join(monthday_climatology, by = "MonthDay")
      working_data$residual <- working_data$target_var - working_data$climatology_mean
      
      # Apply linear detrending
      detrend_data <- working_data %>%
        filter(Year %in% detrend_years) %>%
        filter(!is.na(residual))
      
      detrend_data$time_continuous <- decimal_date(detrend_data$Date)
      
      # Use lm for linear regression
      lm_model <- lm(residual ~ time_continuous, data = detrend_data)
      
      trend_predictions <- predict(lm_model, newdata = detrend_data, se.fit = TRUE)
      detrend_data$trend <- trend_predictions$fit
      detrend_data$trend_se <- trend_predictions$se.fit
      
      p_residual_trend <- ggplot(detrend_data, aes(x = Date)) +
        geom_point(aes(y = residual), alpha = 0.3, size = 0.8, color = "gray50") +
        geom_line(aes(y = trend), color = "blue", size = 1.2) +
        geom_ribbon(aes(ymin = trend - 2 * trend_se, ymax = trend + 2 * trend_se), alpha = 0.2, fill = "blue") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(title = paste(current_var, "Residuals and Linear Trend"),
             subtitle = paste0("Fitted on ", min(detrend_years), "-", max(detrend_years), " data."),
             x = "Date", y = "Residual") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 11))
      
      print(p_residual_trend)
      
      working_data$time_continuous <- decimal_date(working_data$Date)
      working_data$trend <- predict(lm_model, newdata = working_data)
      
      # Calculate the final anomalies
      working_data$anomaly <- working_data$residual - working_data$trend
      
      # Use the standard z-score
      anomaly_mean <- mean(working_data$anomaly, na.rm = TRUE)
      anomaly_sd   <- sd(working_data$anomaly, na.rm = TRUE)
      
      # Calculate the standard z-score
      working_data$z_score <- (working_data$anomaly - anomaly_mean) / anomaly_sd
      
      # Merge the results back into the main data frame
      # Rename columns to include the variable name
      names(working_data)[names(working_data) == "target_var"] <- paste0("actual_", current_var)
      names(working_data)[names(working_data) == "climatology_mean"] <- paste0("climatology_", current_var)
      names(working_data)[names(working_data) == "residual"] <- paste0("residual_", current_var)
      names(working_data)[names(working_data) == "trend"] <- paste0("trend_", current_var)
      names(working_data)[names(working_data) == "anomaly"] <- paste0("anomaly_", current_var)
      names(working_data)[names(working_data) == "z_score"] <- paste0("z_score_", current_var)
      
      # Select only the columns to be merged
      cols_to_merge <- working_data %>% select(Date, starts_with(c("actual", "climatology", "residual", "trend", "anomaly", "z_score")))
      
      # Merge the results using left_join
      all_results_data <- all_results_data %>%
        left_join(cols_to_merge, by = "Date")
      
      # Visualization
      anomaly_col_name <- paste0("anomaly_", current_var)
      z_score_col_name <- paste0("z_score_", current_var)
      
      final_analysis_data <- all_results_data %>%
        filter(Year %in% analysis_years) %>%
        filter(!is.na(!!sym(anomaly_col_name)))
      
      p_anomaly_ts <- ggplot(final_analysis_data, aes(x = Date, y = !!sym(anomaly_col_name))) +
        geom_point(aes(color = abs(!!sym(z_score_col_name)) > extreme_threshold), alpha = 0.5, size = 1) +
        scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                           labels = c("Normal", paste0("Extreme (|z| > ", extreme_threshold, ")")),
                           name = "Anomaly Type") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(title = paste(current_var, "Anomalies (", min(analysis_years), "-", max(analysis_years), ")"),
             x = "Date", y = "Anomaly") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top")
      print(p_anomaly_ts)
      
    }
    
    dev.off()
    
    # Export the merged results
    daily_export_data <- all_results_data %>%
      filter(Year %in% export_years) %>%
      select(Date, Year, Month, Day, MonthDay, DOY, starts_with(c("actual_", "climatology_", "residual_", "anomaly_", "z_score_"))) %>%
      arrange(Date)
    
    output_file_name <- file.path(site_output_dir, paste0("carbon_anomalies_QC_filtered_", file_prefix, ".csv"))
    
    write.csv(daily_export_data, output_file_name, row.names = FALSE)
    cat("  Results saved\n")
    
    writeLines(paste("  SUCCESS: Completed", folder_name), log_con)
    successful <- successful + 1
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    writeLines(paste("  ERROR:", e$message), log_con)
    failed <- failed + 1
  })
  
  # Return to the base directory
  setwd(base_dir)
  
  # Flush the log file
  flush(log_con)
}

# Final summary statistics
cat("\n\n=======================================\n")
cat("Batch processing completed!\n")
cat("Successful:", successful, "\n")
cat("Failed:", failed, "\n")
cat("Skipped:", skipped, "\n")
cat("=======================================\n")

writeLines("\n=======================================", log_con)
writeLines("CARBON ANOMALY BATCH PROCESSING SUMMARY", log_con)
writeLines(paste("Completed at:", Sys.time()), log_con)
writeLines(paste("Successful:", successful), log_con)
writeLines(paste("Failed:", failed), log_con)
writeLines(paste("Skipped:", skipped), log_con)

# Close the log file
close(log_con)

cat("\nSee", log_file, "for details\n")
