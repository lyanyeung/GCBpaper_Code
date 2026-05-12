# ===============================================================================
# PARAMETERIZED HEATWAVE INDEX CALCULATOR (COMBINED ANALYSIS)
# ===============================================================================
# This script calculates heatwave characteristics (frequency, duration, severity)
# for multiple time periods, each with its own reference climatology. It then
# combines the results to construct a single composite heatwave index.
#
# CONFIGURABLE PARAMETERS:
# - Master configuration for shared settings (files, window size, etc.)
# - A list of period-specific configurations for reference/target years.
# ===============================================================================

# Load required libraries
library(copula)
library(VineCopula)
library(kdecopula)
library(scatterplot3d)
library(MASS)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(viridis)
library(patchwork)
library(moments)
library(RColorBrewer)
library(scales)

# ===============================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# ===============================================================================

# ===============================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# ===============================================================================



MASTER_CONFIG <- list(
  # Directory containing the ERA5/ERA5-Land input data for one site.
  # Please replace this with your local data directory.
  working_dir = "~/data/heatwave_inputs/DK_Sor_1996_2020",
  
  # ERA5/ERA5-Land input CSV file for the selected site.
  # Please replace this with the actual input filename in working_dir.
  input_file = "ERA5_or_ERA5Land_input_file_for_this_site.csv",
  
  window_size = 15,
  threshold_quantile = 0.95,
  output_prefix = "~/data/heatwave_outputs/DK-Sor_HW",
  copula_families = c(1:10),
  selection_criterion = "AIC"
)


# List of period-specific configurations to be processed sequentially (UPDATED FOR FOUR PERIODS)
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


# ===============================================================================
# VALIDATION AND SETUP
# ===============================================================================

validate_config <- function(config) {
  cat("===== CONFIGURATION VALIDATION =====\n")
  
  # Validate year ranges
  if (config$reference_end <= config$reference_start) {
    stop("Error: Reference end year must be greater than start year")
  }
  
  if (config$target_end <= config$target_start) {
    stop("Error: Target end year must be greater than start year")
  }
  
  if (config$target_start <= config$reference_end) {
    warning("Warning: Target period overlaps with reference period")
  }
  
  # Check reference period length
  ref_length <- config$reference_end - config$reference_start + 1
  if (ref_length < 20) {
    warning(paste("Warning: Reference period is only", ref_length, "years. Consider using at least 20-30 years."))
  }
  
  cat("Reference period:", config$reference_start, "-", config$reference_end, "(", ref_length, "years)\n")
  cat("Target period:", config$target_start, "-", config$target_end, "(", config$target_end - config$target_start + 1, "years)\n")
  cat("Window size: ±", config$window_size, "days\n")
  cat("Threshold quantile:", config$threshold_quantile, "\n")
  cat("Configuration validated successfully!\n\n")
  
  return(TRUE)
}

# ===============================================================================
# DATA LOADING AND PREPROCESSING
# ===============================================================================

load_and_prepare_data <- function(config) {
  cat("===== DATA LOADING AND PREPROCESSING =====\n")
  
  # Set working directory
  setwd(config$working_dir)
  cat("Working directory:", getwd(), "\n")
  
  # Load data
  cat("Loading data from:", config$input_file, "\n")
  data <- read.csv(config$input_file)
  
  # Date processing
  data$date <- as.Date(data$DATETIME)
  data$year <- year(data$date)
  data$month <- month(data$date)
  data$day <- day(data$date)
  data$doy <- yday(data$date)
  
  cat("Data loaded successfully!\n")
  cat("Date range:", min(data$date, na.rm = TRUE), "to", max(data$date, na.rm = TRUE), "\n")
  cat("Total observations:", nrow(data), "\n\n")
  
  return(data)
}

# ===============================================================================
# THRESHOLD CALCULATION FUNCTIONS
# ===============================================================================

calculate_threshold <- function(target_month, target_day, data_subset, config) {
  # Clean and select relevant data
  data_clean <- data_subset %>%
    dplyr::select(date, year, month, day, doy, temperature_2m_C)
  
  # Create target date for day-of-year calculation
  target_date <- as.Date(paste0("2000-", sprintf("%02d", target_month), "-", sprintf("%02d", target_day)))
  target_doy <- yday(target_date)
  
  # Define window dates
  window_dates <- seq(target_date - config$window_size, target_date + config$window_size, by = "day")
  window_months <- month(window_dates)
  window_days <- day(window_dates)
  
  # Extract data within the window
  window_data <- data_clean %>%
    filter(paste(month, day, sep = "-") %in% paste(window_months, window_days, sep = "-"))
  
  # Calculate threshold
  threshold <- quantile(window_data$temperature_2m_C, probs = config$threshold_quantile, na.rm = TRUE)
  
  return(data.frame(
    month = target_month,
    day = target_day,
    doy = target_doy,
    threshold_95 = as.numeric(threshold),
    n_obs = nrow(window_data)
  ))
}

calculate_daily_thresholds <- function(data, config) {
  cat("===== CALCULATING DAILY THRESHOLDS =====\n")
  cat("Reference period:", config$reference_start, "-", config$reference_end, "\n")
  cat("Window size: ±", config$window_size, "days\n")
  cat("Threshold quantile:", config$threshold_quantile, "\n")
  
  # Filter data for reference period
  data_subset <- data %>%
    filter(year >= config$reference_start & year <= config$reference_end)
  
  cat("Reference data points:", nrow(data_subset), "\n")
  
  # Create all dates for a year (using 2000 as template)
  all_dates <- data.frame(
    date = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), by = "day")
  ) %>%
    mutate(
      month = month(date),
      day = day(date)
    )
  
  # Calculate thresholds for each day of the year
  cat("Calculating thresholds for", nrow(all_dates), "days...\n")
  
  thresholds <- map2_dfr(all_dates$month, all_dates$day,
                         ~calculate_threshold(.x, .y, data_subset, config))
  
  thresholds$date <- as.Date(paste0("2000-",
                                    sprintf("%02d", thresholds$month), "-",
                                    sprintf("%02d", thresholds$day)))
  
  cat("Thresholds calculated successfully!\n")
  cat("Mean threshold:", round(mean(thresholds$threshold_95, na.rm = TRUE), 2), "°C\n")
  cat("Range:", round(min(thresholds$threshold_95, na.rm = TRUE), 2), "to",
      round(max(thresholds$threshold_95, na.rm = TRUE), 2), "°C\n\n")
  
  return(thresholds)
}

# ===============================================================================
# HEATWAVE CHARACTERISTIC CALCULATION FUNCTIONS
# ===============================================================================

calculate_heatwave_frequency <- function(target_year, data, thresholds, config) {
  # Extract extended data (target year ± 1 year to avoid edge effects)
  data_extended <- data %>%
    filter(year >= (target_year - 1) & year <= (target_year + 1)) %>%
    dplyr::select(date, year, month, day, doy, temperature_2m_C) %>%
    arrange(date)
  
  # Join with thresholds
  data_extended <- data_extended %>%
    left_join(
      thresholds %>% dplyr::select(month, day, threshold_95),
      by = c("month", "day")
    )
  
  # Handle February 29th (use February 28th threshold)
  feb28_threshold <- thresholds %>%
    filter(month == 2, day == 28) %>%
    pull(threshold_95)
  
  data_extended <- data_extended %>%
    mutate(threshold_95 = ifelse(
      is.na(threshold_95) & month == 2 & day == 29,
      feb28_threshold,
      threshold_95
    ))
  
  # Flag threshold exceedances
  data_extended$exceed_threshold <- data_extended$temperature_2m_C > data_extended$threshold_95
  
  # Prepare results
  results <- data.frame()
  target_dates <- seq(as.Date(paste0(target_year, "-01-01")),
                      as.Date(paste0(target_year, "-12-31")),
                      by = "day")
  
  # Calculate for each day using sliding window
  for (i in 1:length(target_dates)) {
    current_date <- target_dates[i]
    window_start <- current_date - config$window_size
    window_end <- current_date + config$window_size
    
    window_data <- data_extended %>%
      filter(date >= window_start & date <= window_end) %>%
      arrange(date)
    
    if (nrow(window_data) > 0 && sum(window_data$exceed_threshold, na.rm = TRUE) > 0) {
      rle_result <- rle(window_data$exceed_threshold)
      heatwave_events <- sum(rle_result$values == TRUE, na.rm = TRUE)
      event_lengths <- rle_result$lengths[rle_result$values == TRUE]
      max_duration <- ifelse(length(event_lengths) > 0, max(event_lengths), 0)
      total_days <- sum(window_data$exceed_threshold, na.rm = TRUE)
    } else {
      heatwave_events <- 0
      max_duration <- 0
      total_days <- 0
    }
    
    results <- rbind(results, data.frame(
      date = current_date,
      year = year(current_date),
      month = month(current_date),
      day = day(current_date),
      doy = yday(current_date),
      heatwave_frequency = heatwave_events,
      total_heatwave_days = total_days,
      max_event_duration = max_duration,
      window_start = window_start,
      window_end = window_end
    ))
  }
  
  return(results)
}

calculate_heatwave_severity <- function(target_year, data, thresholds, config) {
  # Extract extended data
  data_extended <- data %>%
    filter(year >= (target_year - 1) & year <= (target_year + 1)) %>%
    dplyr::select(date, year, month, day, doy, temperature_2m_C) %>%
    arrange(date)
  
  # Join with thresholds
  data_extended <- data_extended %>%
    left_join(
      thresholds %>% dplyr::select(month, day, threshold_95),
      by = c("month", "day")
    )
  
  # Handle February 29th
  feb28_threshold <- thresholds %>%
    filter(month == 2, day == 28) %>%
    pull(threshold_95)
  
  data_extended <- data_extended %>%
    mutate(threshold_95 = ifelse(
      is.na(threshold_95) & month == 2 & day == 29,
      feb28_threshold,
      threshold_95
    ))
  
  # Calculate temperature excess
  data_extended <- data_extended %>%
    mutate(
      exceed_threshold = temperature_2m_C > threshold_95,
      temperature_excess = ifelse(exceed_threshold, temperature_2m_C - threshold_95, 0)
    )
  
  # Prepare results
  results <- data.frame()
  target_dates <- seq(as.Date(paste0(target_year, "-01-01")),
                      as.Date(paste0(target_year, "-12-31")),
                      by = "day")
  
  for (i in 1:length(target_dates)) {
    current_date <- target_dates[i]
    window_start <- current_date - config$window_size
    window_end <- current_date + config$window_size
    
    window_data <- data_extended %>%
      filter(date >= window_start & date <= window_end) %>%
      arrange(date)
    
    if (nrow(window_data) > 0) {
      severity <- sum(window_data$temperature_excess, na.rm = TRUE)
      days_exceeding <- sum(window_data$exceed_threshold, na.rm = TRUE)
      max_excess <- ifelse(days_exceeding > 0,
                           max(window_data$temperature_excess, na.rm = TRUE), 0)
      mean_excess <- ifelse(days_exceeding > 0,
                            mean(window_data$temperature_excess[window_data$exceed_threshold], na.rm = TRUE), 0)
      
      if (days_exceeding > 0) {
        rle_result <- rle(window_data$exceed_threshold)
        heatwave_events <- sum(rle_result$values == TRUE, na.rm = TRUE)
      } else {
        heatwave_events <- 0
      }
    } else {
      severity <- 0
      days_exceeding <- 0
      max_excess <- 0
      mean_excess <- 0
      heatwave_events <- 0
    }
    
    results <- rbind(results, data.frame(
      date = current_date,
      year = year(current_date),
      month = month(current_date),
      day = day(current_date),
      doy = yday(current_date),
      severity = severity,
      days_exceeding = days_exceeding,
      max_excess = max_excess,
      mean_excess = mean_excess,
      heatwave_events = heatwave_events
    ))
  }
  
  return(results)
}

# ===============================================================================
# BATCH PROCESSING FUNCTIONS
# ===============================================================================

process_all_years <- function(data, thresholds, config) {
  cat("===== BATCH PROCESSING ALL CHARACTERISTICS =====\n")
  cat("Target period:", config$target_start, "-", config$target_end, "\n")
  cat("Processing", config$target_end - config$target_start + 1, "years...\n\n")
  
  # Initialize storage
  all_frequency_results <- data.frame()
  all_severity_results <- data.frame()
  annual_summary <- data.frame()
  
  # Process each year
  for (year in config$target_start:config$target_end) {
    cat("Processing year", year, "...")
    
    # Calculate frequency (includes duration)
    freq_results <- calculate_heatwave_frequency(year, data, thresholds, config)
    all_frequency_results <- rbind(all_frequency_results, freq_results)
    
    # Calculate severity
    sev_results <- calculate_heatwave_severity(year, data, thresholds, config)
    all_severity_results <- rbind(all_severity_results, sev_results)
    
    # Annual summary
    year_summary <- data.frame(
      year = year,
      total_days = nrow(freq_results),
      days_with_heatwave = sum(freq_results$heatwave_frequency > 0),
      avg_frequency = round(mean(freq_results$heatwave_frequency), 2),
      max_frequency = max(freq_results$heatwave_frequency),
      max_duration = max(freq_results$max_event_duration),
      avg_severity = round(mean(sev_results$severity), 2),
      max_severity = round(max(sev_results$severity), 2)
    )
    
    annual_summary <- rbind(annual_summary, year_summary)
    cat(" Done!\n")
  }
  
  cat("\n===== ANNUAL SUMMARY (for this period) =====\n")
  print(annual_summary)
  
  return(list(
    frequency = all_frequency_results,
    severity = all_severity_results,
    summary = annual_summary
  ))
}

# ===============================================================================
# COPULA MODELING AND INDEX CONSTRUCTION
# ===============================================================================

construct_heatwave_index <- function(all_results, config) {
  cat("\n===== CONSTRUCTING HEATWAVE INDEX (on combined data) =====\n")
  
  # Merge datasets
  characteristics_data <- all_results$frequency %>%
    dplyr::select(date, year, month, day, doy,
                  heatwave_frequency, max_event_duration,
                  window_start, window_end) %>%
    left_join(
      all_results$severity %>%
        dplyr::select(date, severity, days_exceeding),
      by = "date"
    )
  
  cat("Data merged successfully!\n")
  cat("Total observations:", nrow(characteristics_data), "\n")
  
  # Data preprocessing for copula
  data_heatwave <- characteristics_data %>%
    filter(heatwave_frequency > 0 | max_event_duration > 0 | severity > 0) %>%
    dplyr::select(date, year, month, doy, heatwave_frequency, max_event_duration, severity)
  
  cat("Days with heatwave events:", nrow(data_heatwave), "\n")
  cat("Days without heatwave events:", nrow(characteristics_data) - nrow(data_heatwave), "\n")
  
  if (nrow(data_heatwave) == 0) {
    stop("Error: No heatwave events found in the data!")
  }
  
  # Marginal distribution transformation
  cat("\nTransforming to uniform margins...\n")
  u_freq <- rank(data_heatwave$heatwave_frequency) / (nrow(data_heatwave) + 1)
  u_duration <- rank(data_heatwave$max_event_duration) / (nrow(data_heatwave) + 1)
  u_severity <- rank(data_heatwave$severity) / (nrow(data_heatwave) + 1)
  
  u_data <- data.frame(
    u_freq = u_freq,
    u_duration = u_duration,
    u_severity = u_severity
  )
  
  # Vine copula modeling
  cat("Fitting vine copula model...\n")
  vine_fit <- RVineStructureSelect(u_data,
                                   familyset = config$copula_families,
                                   selectioncrit = config$selection_criterion,
                                   indeptest = TRUE)
  
  cat("Vine copula fitted successfully!\n")
  
  # Calculate joint CDF and transform
  cat("Calculating joint CDF...\n")
  joint_cdf <- RVineCDF(u_data, vine_fit)
  
  # Transform to normal scale
  normal_transformed <- qnorm(pmin(pmax(joint_cdf, 1e-6), 1 - 1e-6))
  normal_transformed <- scale(normal_transformed)[,1]
  
  # Construct complete index
  data_all <- characteristics_data %>%
    dplyr::select(date, year, month, doy, heatwave_frequency, max_event_duration, severity, days_exceeding, window_start, window_end) %>%
    mutate(has_heatwave = ifelse(heatwave_frequency > 0 | max_event_duration > 0 | severity > 0, 1, 0))
  
  n_total <- nrow(data_all)
  heatwave_indices <- which(data_all$has_heatwave == 1)
  no_heatwave_indices <- which(data_all$has_heatwave == 0)
  
  normal_index <- numeric(n_total)
  normal_index[heatwave_indices] <- normal_transformed
  normal_index[no_heatwave_indices] <- -4 # Assign a low value for non-heatwave days
  
  # Handle extreme values that might fall below the non-heatwave threshold
  problem_rows <- which(data_all$has_heatwave == 1 & normal_index < -4)
  if (length(problem_rows) > 0) {
    cat("Adjusting", length(problem_rows), "extreme values...\n")
    normal_index[problem_rows] <- -4
  }
  
  # Final dataset
  final_data <- data_all %>%
    mutate(normal_index = normal_index)
  
  cat("Heatwave index constructed successfully!\n")
  cat("Index range:", round(min(final_data$normal_index), 2), "to", round(max(final_data$normal_index), 2), "\n")
  
  return(list(
    data = final_data,
    vine_model = vine_fit,
    summary_stats = list(
      total_days = nrow(final_data),
      heatwave_days = sum(final_data$has_heatwave),
      index_range = range(final_data$normal_index)
    )
  ))
}

# ===============================================================================
# OUTPUT AND VISUALIZATION FUNCTIONS
# ===============================================================================

save_results <- function(final_results, final_config, period_configs) {
  cat("\n===== SAVING RESULTS =====\n")
  
  # Generate output filename based on the full analysis range
  output_filename <- paste0(final_config$output_prefix, "_", final_config$target_start, "_", final_config$target_end, ".csv")
  
  # Save main results
  write.csv(final_results$data, file = output_filename, row.names = FALSE)
  cat("Results saved to:", output_filename, "\n")
  
  # Save summary statistics
  summary_filename <- paste0(final_config$output_prefix, "_summary_", final_config$target_start, "_", final_config$target_end, ".txt")
  
  sink(summary_filename)
  cat("HEATWAVE INDEX CALCULATION SUMMARY\n")
  cat("=====================================\n\n")
  cat("Overall Target Period:", final_config$target_start, "-", final_config$target_end, "\n\n")
  cat("Configuration Details:\n")
  
  # Detail the reference period used for each target period
  for (i in 1:length(period_configs)) {
    p_conf <- period_configs[[i]]
    cat(sprintf("- For Target Years %d-%d, Reference Period %d-%d was used.\n", 
                p_conf$target_start, p_conf$target_end, p_conf$reference_start, p_conf$reference_end))
  }
  
  cat("\n- Window size: ±", final_config$window_size, "days\n")
  cat("- Threshold quantile:", final_config$threshold_quantile, "\n\n")
  
  cat("Results:\n")
  cat("- Total days analyzed:", final_results$summary_stats$total_days, "\n")
  cat("- Days with heatwave events:", final_results$summary_stats$heatwave_days, "\n")
  cat("- Index range:", paste(round(final_results$summary_stats$index_range, 2), collapse = " to "), "\n")
  sink()
  
  cat("Summary saved to:", summary_filename, "\n")
}

# ===============================================================================
# MAIN EXECUTION FUNCTION
# ===============================================================================

run_combined_heatwave_analysis <- function(master_config, period_configs) {
  cat("===============================================================================\n")
  cat("COMBINED HEATWAVE INDEX ANALYSIS\n")
  cat("===============================================================================\n\n")
  
  # Step 1: Load and prepare data once
  data <- load_and_prepare_data(master_config)
  
  # Initialize containers for results from all periods
  combined_frequency_results <- data.frame()
  combined_severity_results <- data.frame()
  
  # Step 2: Loop through each period configuration
  for (i in 1:length(period_configs)) {
    period_name <- names(period_configs)[i]
    period_conf <- period_configs[[i]]
    
    cat(sprintf("\n--- Processing %s (Target: %d-%d, Reference: %d-%d) ---\n\n",
                period_name, period_conf$target_start, period_conf$target_end,
                period_conf$reference_start, period_conf$reference_end))
    
    # Merge master and period-specific configs
    current_config <- c(master_config, period_conf)
    
    # Validate the current configuration
    validate_config(current_config)
    
    # Calculate thresholds for the current reference period
    thresholds <- calculate_daily_thresholds(data, current_config)
    
    # Process all years in the current target period
    period_results <- process_all_years(data, thresholds, current_config)
    
    # Append results to the combined data frames
    combined_frequency_results <- rbind(combined_frequency_results, period_results$frequency)
    combined_severity_results <- rbind(combined_severity_results, period_results$severity)
  }
  
  # Step 3: Combine all results for final index construction
  combined_all_results <- list(
    frequency = combined_frequency_results,
    severity = combined_severity_results
  )
  
  # Create a final configuration object for the index and saving steps
  # This represents the entire analysis span
  final_config <- master_config
  final_config$target_start <- min(sapply(period_configs, `[[`, "target_start"))
  final_config$target_end <- max(sapply(period_configs, `[[`, "target_end"))
  
  # Step 4: Construct heatwave index using the combined data
  final_results <- construct_heatwave_index(combined_all_results, final_config)
  
  # Step 5: Save the final, combined results
  save_results(final_results, final_config, period_configs)
  
  cat("\n===============================================================================\n")
  cat("COMBINED ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("===============================================================================\n")
  
  return(final_results)
}

# ===============================================================================
# RUN THE ANALYSIS
# ===============================================================================

# Execute the combined analysis with the configurations defined at the top
results <- run_combined_heatwave_analysis(MASTER_CONFIG, PERIOD_CONFIGS)
