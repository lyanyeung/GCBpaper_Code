# =============================================================================
# Carbon Uptake Interval Analysis
#
# Purpose:
#   This script calculates site-specific climatological seasonal cycles for
#   NEE, GPP, and ecosystem respiration, and identifies carbon uptake intervals
#   as contiguous periods when climatological NEE < 0.
#
# Method:
#   Climatology is calculated using a month-day based 31-day moving window.
#   This avoids relying directly on calendar-year DOY and provides consistent
#   handling of leap years.
#
# Input data structure:
#   data/
#     Folder_List_All.csv
#     SITE_YYYY_YYYY/
#       SITE.csv
#
# Example:
#   data/
#     Folder_List_All.csv
#     AU_How_2001_2014/
#       AU_How.csv
#
# Required columns in each site CSV:
#   TIMESTAMP
#   NEE_VUT_REF
#   GPP_NT_VUT_REF
#   RECO_NT_VUT_REF
#
# Main output:
#   All_Sites_Carbon_Uptake_Intervals.csv
#
# Notes:
#   - Negative NEE indicates net carbon uptake.
#   - The baseline period for each site is extracted from the folder name,
#     e.g. AU_How_2001_2014 uses 2001-2014 as the baseline years.
# =============================================================================


# =============================================================================
# Section 0: User Configuration
# =============================================================================

BATCH_CONFIG <- list(
  
  # ---------------------------------------------------------------------------
  # Paths
  # ---------------------------------------------------------------------------
  # Set this to the folder containing Folder_List_All.csv and all site folders.
  # For local use, replace "data" with your own data directory, for example:
  # base_dir = "E:/ECTower_15_0.7"
  
  base_dir = "data",
  
  # CSV file listing all site folders.
  # The file should contain a column named "folder".
  folder_list_file = "Folder_List_All.csv",
  
  # Output directory. It will be created automatically if it does not exist.
  output_base_dir = "outputs/carbon_uptake_intervals",
  
  # Log file name saved inside output_base_dir.
  log_file = "batch_climatology_log.txt",
  
  # If TRUE, individual site-level climatology PDFs are deleted after they are
  # combined into a single PDF.
  delete_individual_pdfs = TRUE,
  
  
  # ---------------------------------------------------------------------------
  # Analysis settings
  # ---------------------------------------------------------------------------
  # Variables for which climatology will be calculated.
  target_vars = c("NEE_VUT_REF", "GPP_NT_VUT_REF", "RECO_NT_VUT_REF"),
  
  # Variable used to define carbon uptake intervals.
  # By convention, NEE < 0 indicates net ecosystem carbon uptake.
  uptake_variable = "NEE_VUT_REF",
  
  # Moving-window size for climatology calculation.
  # This must be an odd number so that the window is centred on each calendar day.
  window_size = 31
)

# =============================================================================
# Section 1: Library Loading and Setup
# =============================================================================
cat("========== Batch Climatology & Uptake Intervals Analysis (Month-day workflow) ==========\n")

library(tidyverse)
library(lubridate)
library(ggplot2)
library(qpdf)

if (BATCH_CONFIG$window_size %% 2 == 0) stop("Error: window_size must be an odd number.")

# =============================================================================
# Section 2: Helper Functions
# =============================================================================

extract_site_info <- function(folder_name) {
  clean_name <- sub("_\\d{4}_\\d{4}$", "", folder_name)
  country_code <- substr(clean_name, 1, 2)
  underscore_pos <- regexpr("_", clean_name)[1]
  site_part <- if (underscore_pos > 0) substr(clean_name, underscore_pos + 1, nchar(clean_name)) else substr(clean_name, 3, nchar(clean_name))
  full_code <- paste0(country_code, "-", site_part)
  return(list(full_code = full_code, clean_name = clean_name))
}

# Month-day based mean calculation
calculate_monthday_stats <- function(month, day, data_for_calc, window_size) {
  window_offset <- (window_size - 1) / 2
  target_monthday <- sprintf("%02d-%02d", month, day)
  
  # Get all dates in the window
  window_dates <- c()
  
  for (offset in -window_offset:window_offset) {
    # Use 2020 (a leap year) as reference for date arithmetic
    window_date <- as.Date(paste(2020, month, day, sep="-")) + offset
    window_month <- month(window_date)
    window_day <- day(window_date)
    window_monthday <- sprintf("%02d-%02d", window_month, window_day)
    
    # Special handling for Feb 29
    if (window_monthday == "02-29") {
      # Only include if we have Feb 29 data
      if (any(data_for_calc$MonthDay == "02-29")) {
        window_dates <- c(window_dates, window_monthday)
      }
    } else {
      window_dates <- c(window_dates, window_monthday)
    }
  }
  
  window_subset <- data_for_calc %>% filter(MonthDay %in% window_dates)
  
  if (nrow(window_subset) > 0) {
    return(data.frame(
      mean = mean(window_subset$value, na.rm = TRUE),
      sd = sd(window_subset$value, na.rm = TRUE)
    ))
  } else {
    return(data.frame(mean = NA, sd = NA))
  }
}

# Process a single site
process_single_site <- function(folder_name, site_baseline_years, config, log_con) {
  
  site_info <- extract_site_info(folder_name)
  cat("\n========================================\n")
  cat("Processing Site:", site_info$full_code, "\n")
  cat("Using Baseline:", min(site_baseline_years), "-", max(site_baseline_years), "\n")
  writeLines(paste("\n=== Processing", folder_name, "==="), log_con)
  
  input_filename <- paste0(site_info$clean_name, ".csv")
  file_path <- file.path(config$base_dir, folder_name, input_filename)
  
  tryCatch({
    # --- Step 1: Read and Prepare Data ---
    if (!file.exists(file_path)) stop(paste("File not found:", file_path))
    data_raw <- read.csv(file_path, stringsAsFactors = FALSE, na.strings = "-9999")
    
    data <- data_raw %>%
      select(TIMESTAMP, all_of(config$target_vars)) %>%
      mutate(Date = as.Date(as.character(TIMESTAMP), format = "%Y%m%d")) %>%
      filter(!is.na(Date)) %>%
      mutate(
        Year = year(Date),
        Month = month(Date),
        Day = day(Date),
        DOY = yday(Date),  # Keep for plotting
        MonthDay = sprintf("%02d-%02d", Month, Day)
      )
    
    # --- Step 2: Calculate Climatology using Month-day ---
    all_climatologies <- list()
    
    for (current_var in config$target_vars) {
      var_data <- data %>%
        select(Year, Month, Day, MonthDay, value = all_of(current_var)) %>%
        filter(!is.na(value), Year %in% site_baseline_years)
      
      # Get all unique month-days in the data
      unique_monthdays <- unique(var_data$MonthDay)
      climatology_for_var <- data.frame()
      
      for (md in unique_monthdays) {
        month_val <- as.integer(substr(md, 1, 2))
        day_val <- as.integer(substr(md, 4, 5))
        stats <- calculate_monthday_stats(month_val, day_val, var_data, config$window_size)
        stats$MonthDay <- md
        stats$Month <- month_val
        stats$Day <- day_val
        climatology_for_var <- rbind(climatology_for_var, stats)
      }
      
      # Convert to DOY for plotting (using a reference year)
      climatology_for_var$DOY_plotting <- yday(as.Date(paste("2020", 
                                                             climatology_for_var$Month, 
                                                             climatology_for_var$Day, 
                                                             sep="-")))
      climatology_for_var$Variable <- current_var
      all_climatologies[[current_var]] <- climatology_for_var
    }
    
    final_climatology_data <- bind_rows(all_climatologies)
    
    # --- Step 3: Create and Save Climatology Plot ---
    climatology_plot <- ggplot(final_climatology_data, aes(x = DOY_plotting, y = mean)) +
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, fill = "steelblue") +
      geom_line(color = "navy", linewidth = 1) +
      facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
      labs(title = paste("Mean Annual Cycle (Climatology) for Site:", site_info$full_code),
           subtitle = paste("Baseline:", min(site_baseline_years), "-", max(site_baseline_years), 
                            "| Window:", config$window_size, "days | Month-day method"),
           x = "Day of Year", y = "Daily Mean Value") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            strip.text = element_text(face = "bold", size = 12))
    
    output_pdf_path <- file.path(config$output_base_dir, paste0(site_info$full_code, "_climatology.pdf"))
    ggsave(output_pdf_path, climatology_plot, width = 8, height = 10, device = "pdf")
    cat("  - Site climatology plot saved to:", output_pdf_path, "\n")
    
    # --- Step 4: Calculate Carbon Uptake Intervals ---
    nee_climatology <- final_climatology_data %>% 
      filter(Variable == config$uptake_variable) %>%
      arrange(DOY_plotting)
    
    # Create a complete year sequence to handle gaps
    all_days <- data.frame(DOY_plotting = 1:366)
    nee_complete <- all_days %>%
      left_join(nee_climatology %>% select(DOY_plotting, mean, MonthDay), by = "DOY_plotting") %>%
      filter(!is.na(mean))  # Remove days without data (e.g., Feb 29 in non-leap years)
    
    # Identify contiguous negative blocks
    nee_complete <- nee_complete %>%
      mutate(
        is_negative = mean < 0,
        block_id = cumsum(is_negative != lag(is_negative, default = first(is_negative)))
      )
    
    # Find intervals
    intervals <- nee_complete %>%
      filter(is_negative == TRUE) %>%
      group_by(block_id) %>%
      summarise(
        start_doy = min(DOY_plotting), 
        end_doy = max(DOY_plotting),
        start_monthday = first(MonthDay),
        end_monthday = last(MonthDay),
        .groups = 'drop'
      ) %>%
      select(start_doy, end_doy, start_monthday, end_monthday)
    
    # Handle wrap-around seasons
    if (nrow(intervals) > 1) {
      # Check if we have a season that wraps around the year
      first_interval <- intervals[1, ]
      last_interval <- intervals[nrow(intervals), ]
      
      # If the first interval starts at DOY 1 and the last ends at max DOY
      max_doy <- max(nee_complete$DOY_plotting)
      if (first_interval$start_doy == 1 && last_interval$end_doy == max_doy) {
        cat("  - Found a wrap-around season. Merging intervals.\n")
        
        # Create merged interval
        merged_interval <- data.frame(
          start_doy = last_interval$start_doy,
          end_doy = first_interval$end_doy,
          start_monthday = last_interval$start_monthday,
          end_monthday = first_interval$end_monthday
        )
        
        # Remove first and last intervals, add merged
        intervals <- intervals[-c(1, nrow(intervals)), ]
        intervals <- rbind(intervals, merged_interval)
      }
    }
    
    if (nrow(intervals) == 0) {
      cat("  - No carbon uptake period found for this site.\n")
      writeLines(paste("  No uptake intervals found for", folder_name), log_con)
      return(list(success = TRUE, site_code = site_info$full_code, intervals = NULL))
    }
    
    # Add site code
    intervals$site_code <- site_info$full_code
    cat("  - Found", nrow(intervals), "uptake interval(s).\n")
    
    writeLines(paste("  SUCCESS: Completed", folder_name, "with", nrow(intervals), "intervals"), log_con)
    return(list(success = TRUE, site_code = site_info$full_code, intervals = intervals))
    
  }, error = function(e) {
    cat("  - ERROR:", e$message, "\n")
    writeLines(paste("  ERROR for", folder_name, ":", e$message), log_con)
    return(list(success = FALSE, site_code = site_info$full_code, intervals = NULL))
  })
}

# =============================================================================
# Section 3: Main Execution
# =============================================================================

if (!dir.exists(BATCH_CONFIG$output_base_dir)) {
  dir.create(BATCH_CONFIG$output_base_dir, recursive = TRUE)
}

folder_list_path <- file.path(BATCH_CONFIG$base_dir, BATCH_CONFIG$folder_list_file)
if (!file.exists(folder_list_path)) stop("Folder list file not found at:", folder_list_path)

folder_list <- read.csv(folder_list_path)
log_file_path <- file.path(BATCH_CONFIG$output_base_dir, BATCH_CONFIG$log_file)
log_con <- file(log_file_path, open = "wt")
writeLines(paste("Batch processing started at:", Sys.time()), log_con)
writeLines(paste("Using Month-day workflow for leap year handling"), log_con)

# Initialize storage for results
all_site_intervals <- list()
successful <- 0
failed <- 0

# Process each site
for (i in 1:nrow(folder_list)) {
  folder_name <- folder_list$folder[i]
  year_matches <- stringr::str_match(folder_name, "_(\\d{4})_(\\d{4})$")
  
  if (is.na(year_matches[1,1])) {
    cat("\nWARNING: Could not extract years from '", folder_name, "'. Skipping.\n", sep="")
    writeLines(paste("WARNING: Skipping", folder_name, "- cannot extract years"), log_con)
    failed <- failed + 1
    next
  }
  
  site_baseline_years <- as.numeric(year_matches[1, 2]):as.numeric(year_matches[1, 3])
  
  # Process the site
  result <- process_single_site(folder_name, site_baseline_years, BATCH_CONFIG, log_con)
  
  if (result$success) {
    successful <- successful + 1
    if (!is.null(result$intervals)) {
      all_site_intervals[[folder_name]] <- result$intervals
    }
  } else {
    failed <- failed + 1
  }
  
  flush(log_con)
}

# =============================================================================
# Section 4: Export Carbon Uptake Intervals CSV
# =============================================================================
cat("\n###############################################\n")
cat("Generating Carbon Uptake Intervals CSV file...\n")

if (length(all_site_intervals) > 0) {
  # Combine all intervals into one dataframe
  uptake_intervals_df <- bind_rows(all_site_intervals) %>%
    select(site_code, start_doy, end_doy, start_monthday, end_monthday) %>%
    arrange(site_code, start_doy)
  
  output_csv_path <- file.path(BATCH_CONFIG$output_base_dir, "All_Sites_Carbon_Uptake_Intervals.csv")
  write.csv(uptake_intervals_df, output_csv_path, row.names = FALSE)
  cat("SUCCESS: Uptake intervals summary saved to:", output_csv_path, "\n")
  writeLines(paste("Uptake intervals CSV saved:", output_csv_path), log_con)
} else {
  cat("No uptake intervals were found for any site.\n")
  writeLines("No uptake intervals found for any site", log_con)
}

# =============================================================================
# Section 5: Combine All PDFs
# =============================================================================
cat("\n###############################################\n")
cat("Combining all site PDFs...\n")

individual_pdf_files <- list.files(
  path = BATCH_CONFIG$output_base_dir,
  pattern = "_climatology\\.pdf$",
  full.names = TRUE
)

if (length(individual_pdf_files) > 0) {
  combined_pdf_path <- file.path(BATCH_CONFIG$output_base_dir, "All_Sites_Combined_Climatology.pdf")
  tryCatch({
    qpdf::pdf_combine(input = individual_pdf_files, output = combined_pdf_path)
    cat("SUCCESS: Combined", length(individual_pdf_files), "plots into:", combined_pdf_path, "\n")
    writeLines(paste("Combined PDFs into:", combined_pdf_path), log_con)
    
    if (BATCH_CONFIG$delete_individual_pdfs) {
      file.remove(individual_pdf_files)
      cat("Cleaned up individual site PDF files.\n")
      writeLines("Deleted individual PDF files", log_con)
    }
  }, error = function(e) {
    cat("ERROR during final PDF combination:", e$message, "\n")
    writeLines(paste("ERROR combining PDFs:", e$message), log_con)
  })
} else {
  cat("No individual site PDFs found to combine.\n")
  writeLines("No PDFs found to combine", log_con)
}

# --- Final Summary ---
summary_message <- paste0(
  "\n\nBATCH PROCESSING COMPLETE\n",
  "=====================\n",
  "Method: Month-day workflow (preserves Feb 29)\n",
  "Successful: ", successful, "\n",
  "Failed: ", failed, "\n",
  "Results saved in: ", BATCH_CONFIG$output_base_dir, "\n"
)
cat(summary_message)
writeLines(summary_message, log_con)
close(log_con)

cat("\nCheck", BATCH_CONFIG$log_file, "for detailed information.\n")
