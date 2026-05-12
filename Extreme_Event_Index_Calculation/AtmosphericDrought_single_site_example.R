# ===============================================================================
# PARAMETERIZED VPD DROUGHT INDEX CALCULATOR (COMBINED ANALYSIS)
# ===============================================================================
# This script calculates VPD drought characteristics for multiple time periods,
# each with its own reference climatology. It then combines the results
# to construct a single composite drought index using vine copulas.
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
library(purrr)

# ===============================================================================
# CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
# ===============================================================================

# Master configuration for settings common to all analyses
MASTER_CONFIG <- list(
  # Working directory
  working_dir = "C:/Users/ly3n24/OneDrive - University of Southampton/Fluxdata/ECTower/DK_Sor_1996_2020",
  
  # Input data file
  input_file = "DK-Sor_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv",
  
  # Moving window size (days on each side of target day)
  window_size = 15,
  
  # Quantile for threshold calculation (95th percentile for high VPD)
  threshold_quantile = 0.95,
  
  # Output file prefix
  output_prefix = "DK-Sor_VPD_Drought",
  
  # Vine copula family set
  copula_families = c(1:10),
  
  # Selection criterion for vine copula
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

validate_vpd_config <- function(config) {
  cat("===== VPD DROUGHT CONFIGURATION VALIDATION =====\n")
  
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
  cat("VPD threshold quantile:", config$threshold_quantile, "\n")
  cat("Configuration validated successfully!\n\n")
  
  return(TRUE)
}

# ===============================================================================
# DATA LOADING AND PREPROCESSING
# ===============================================================================

load_and_prepare_vpd_data <- function(config) {
  cat("===== VPD DATA LOADING AND PREPROCESSING =====\n")
  
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
# VPD THRESHOLD CALCULATION FUNCTIONS
# ===============================================================================

calculate_vpd_threshold <- function(target_month, target_day, data_subset, config) {
  # Select relevant data
  data_clean <- data_subset %>%
    dplyr::select(date, year, month, day, doy, VPD_kPa)
  
  # Create a template date for window calculation
  target_date <- as.Date(paste0("2000-", sprintf("%02d", target_month), "-", sprintf("%02d", target_day)))
  target_doy <- yday(target_date)
  
  # Define window
  window_dates <- seq(target_date - config$window_size, target_date + config$window_size, by = "day")
  window_months <- month(window_dates)
  window_days <- day(window_dates)
  
  # Filter data within the window across all years in the subset
  window_data <- data_clean %>%
    filter(paste(month, day, sep = "-") %in% paste(window_months, window_days, sep = "-"))
  
  # Calculate VPD threshold
  vpd_threshold <- quantile(window_data$VPD_kPa, probs = config$threshold_quantile, na.rm = TRUE)
  
  return(data.frame(
    month = target_month,
    day = target_day,
    doy = target_doy,
    vpd_threshold_95 = as.numeric(vpd_threshold),
    n_obs = nrow(window_data)
  ))
}

calculate_daily_vpd_thresholds <- function(data, config) {
  cat("===== CALCULATING DAILY VPD THRESHOLDS =====\n")
  cat("Reference period:", config$reference_start, "-", config$reference_end, "\n")
  cat("Window size: ±", config$window_size, "days\n")
  cat("VPD threshold quantile:", config$threshold_quantile, "\n")
  
  # Filter data for the specified reference period
  data_subset <- data %>%
    filter(year >= config$reference_start & year <= config$reference_end)
  
  cat("Reference data points:", nrow(data_subset), "\n")
  
  # Create a template for all days of a year
  all_dates <- data.frame(date = seq(as.Date("2000-01-01"), as.Date("2000-12-31"), by = "day")) %>%
    mutate(month = month(date), day = day(date))
  
  # Calculate thresholds for each day
  cat("Calculating VPD thresholds for", nrow(all_dates), "days...\n")
  vpd_thresholds <- map2_dfr(all_dates$month, all_dates$day,
                             ~calculate_vpd_threshold(.x, .y, data_subset, config))
  
  cat("VPD thresholds calculated successfully!\n")
  cat("Mean VPD threshold:", round(mean(vpd_thresholds$vpd_threshold_95, na.rm = TRUE), 3), "kPa\n")
  cat("Range:", round(min(vpd_thresholds$vpd_threshold_95, na.rm = TRUE), 3), "to",
      round(max(vpd_thresholds$vpd_threshold_95, na.rm = TRUE), 3), "kPa\n\n")
  
  return(vpd_thresholds)
}

# ===============================================================================
# VPD DROUGHT CHARACTERISTIC CALCULATION FUNCTIONS
# ===============================================================================

calculate_drought_frequency <- function(target_year, data, vpd_thresholds, config) {
  # Widen data to include adjacent years to handle window edge effects
  data_extended <- data %>%
    filter(year >= (target_year - 1) & year <= (target_year + 1)) %>%
    dplyr::select(date, year, month, day, doy, VPD_kPa) %>%
    arrange(date)
  
  # Join with daily VPD thresholds
  data_extended <- data_extended %>%
    left_join(vpd_thresholds %>% dplyr::select(month, day, vpd_threshold_95), by = c("month", "day"))
  
  # Handle leap years by using Feb 28th threshold for Feb 29th
  feb28_vpd_threshold <- vpd_thresholds %>% filter(month == 2, day == 28) %>% pull(vpd_threshold_95)
  data_extended <- data_extended %>%
    mutate(vpd_threshold_95 = ifelse(is.na(vpd_threshold_95) & month == 2 & day == 29, feb28_vpd_threshold, vpd_threshold_95))
  
  # A day is in drought if its VPD exceeds the daily threshold
  data_extended$exceed_threshold <- data_extended$VPD_kPa > data_extended$vpd_threshold_95
  
  # Loop through each day of the target year to calculate characteristics within its moving window
  results <- data.frame()
  target_dates <- seq(as.Date(paste0(target_year, "-01-01")), as.Date(paste0(target_year, "-12-31")), by = "day")
  
  for (i in 1:length(target_dates)) {
    current_date <- target_dates[i]
    window_start <- current_date - config$window_size
    window_end <- current_date + config$window_size
    
    window_data <- data_extended %>% filter(date >= window_start & date <= window_end)
    
    if (nrow(window_data) > 0 && sum(window_data$exceed_threshold, na.rm = TRUE) > 0) {
      rle_result <- rle(window_data$exceed_threshold)
      drought_events <- sum(rle_result$values == TRUE, na.rm = TRUE)
      event_lengths <- rle_result$lengths[rle_result$values == TRUE]
      max_duration <- ifelse(length(event_lengths) > 0, max(event_lengths), 0)
      total_days <- sum(window_data$exceed_threshold, na.rm = TRUE)
    } else {
      drought_events <- 0
      max_duration <- 0
      total_days <- 0
    }
    
    results <- rbind(results, data.frame(
      date = current_date,
      year = year(current_date),
      month = month(current_date),
      day = day(current_date),
      doy = yday(current_date),
      drought_frequency = drought_events,
      total_drought_days = total_days,
      max_event_duration = max_duration,
      window_start = window_start,
      window_end = window_end
    ))
  }
  return(results)
}

calculate_drought_severity <- function(target_year, data, vpd_thresholds, config) {
  data_extended <- data %>%
    filter(year >= (target_year - 1) & year <= (target_year + 1)) %>%
    dplyr::select(date, year, month, day, doy, VPD_kPa) %>%
    arrange(date)
  
  data_extended <- data_extended %>%
    left_join(vpd_thresholds %>% dplyr::select(month, day, vpd_threshold_95), by = c("month", "day"))
  
  feb28_vpd_threshold <- vpd_thresholds %>% filter(month == 2, day == 28) %>% pull(vpd_threshold_95)
  data_extended <- data_extended %>%
    mutate(vpd_threshold_95 = ifelse(is.na(vpd_threshold_95) & month == 2 & day == 29, feb28_vpd_threshold, vpd_threshold_95))
  
  # Calculate severity as the cumulative VPD amount exceeding the threshold
  data_extended <- data_extended %>%
    mutate(
      exceed_threshold = VPD_kPa > vpd_threshold_95,
      vpd_excess = ifelse(exceed_threshold, VPD_kPa - vpd_threshold_95, 0)
    )
  
  results <- data.frame()
  target_dates <- seq(as.Date(paste0(target_year, "-01-01")), as.Date(paste0(target_year, "-12-31")), by = "day")
  
  for (i in 1:length(target_dates)) {
    current_date <- target_dates[i]
    window_start <- current_date - config$window_size
    window_end <- current_date + config$window_size
    
    window_data <- data_extended %>% filter(date >= window_start & date <= window_end)
    
    if (nrow(window_data) > 0) {
      severity <- sum(window_data$vpd_excess, na.rm = TRUE)
      days_exceeding <- sum(window_data$exceed_threshold, na.rm = TRUE)
      max_excess <- ifelse(days_exceeding > 0, max(window_data$vpd_excess, na.rm = TRUE), 0)
      mean_excess <- ifelse(days_exceeding > 0, mean(window_data$vpd_excess[window_data$exceed_threshold], na.rm = TRUE), 0)
    } else {
      severity <- 0; days_exceeding <- 0; max_excess <- 0; mean_excess <- 0
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
      mean_excess = mean_excess
    ))
  }
  return(results)
}

# ===============================================================================
# BATCH PROCESSING FUNCTIONS
# ===============================================================================

process_all_drought_years <- function(data, vpd_thresholds, config) {
  cat("===== BATCH PROCESSING VPD DROUGHT CHARACTERISTICS =====\n")
  cat("Target period:", config$target_start, "-", config$target_end, "\n")
  
  all_frequency_results <- data.frame()
  all_severity_results <- data.frame()
  
  for (year in config$target_start:config$target_end) {
    cat("Processing year", year, "...\n")
    freq_results <- calculate_drought_frequency(year, data, vpd_thresholds, config)
    sev_results <- calculate_drought_severity(year, data, vpd_thresholds, config)
    all_frequency_results <- rbind(all_frequency_results, freq_results)
    all_severity_results <- rbind(all_severity_results, sev_results)
  }
  
  cat("\nBatch processing for period complete.\n")
  return(list(frequency = all_frequency_results, severity = all_severity_results))
}

# ===============================================================================
# COPULA MODELING AND DROUGHT INDEX CONSTRUCTION
# ===============================================================================

construct_drought_index <- function(all_results, config) {
  cat("\n===== CONSTRUCTING VPD DROUGHT INDEX (on combined data) =====\n")
  
  # Merge frequency and severity results
  characteristics_data <- all_results$frequency %>%
    left_join(all_results$severity %>% dplyr::select(date, severity, days_exceeding), by = "date")
  
  cat("Data merged successfully! Total observations:", nrow(characteristics_data), "\n")
  
  # Filter for days that are part of a drought event
  data_drought <- characteristics_data %>%
    filter(drought_frequency > 0 | max_event_duration > 0 | severity > 0) %>%
    dplyr::select(date, drought_frequency, max_event_duration, severity)
  
  if (nrow(data_drought) == 0) {
    stop("Error: No drought events found in the combined data!")
  }
  cat("Days with drought events:", nrow(data_drought), "\n")
  
  # Transform marginal distributions to uniform scale [0,1]
  cat("Transforming to uniform margins...\n")
  u_data <- data.frame(
    u_freq = rank(data_drought$drought_frequency) / (nrow(data_drought) + 1),
    u_duration = rank(data_drought$max_event_duration) / (nrow(data_drought) + 1),
    u_severity = rank(data_drought$severity) / (nrow(data_drought) + 1)
  )
  
  # Fit the vine copula model
  cat("Fitting vine copula model...\n")
  vine_fit_drought <- RVineStructureSelect(u_data,
                                           familyset = config$copula_families,
                                           selectioncrit = config$selection_criterion,
                                           indeptest = TRUE)
  cat("Vine copula fitted successfully!\n")
  
  # Calculate the joint probability (CDF)
  cat("Calculating joint CDF...\n")
  joint_cdf <- RVineCDF(u_data, vine_fit_drought)
  
  # Transform the CDF values to a standard normal scale for the index
  normal_transformed <- qnorm(pmin(pmax(joint_cdf, 1e-6), 1 - 1e-6))
  normal_transformed <- scale(normal_transformed)[, 1]
  
  # Integrate the index back into the full dataset
  data_all <- characteristics_data %>%
    mutate(has_drought = ifelse(drought_frequency > 0 | max_event_duration > 0 | severity > 0, 1, 0))
  
  normal_index <- numeric(nrow(data_all))
  normal_index[data_all$has_drought == 1] <- normal_transformed
  normal_index[data_all$has_drought == 0] <- -4 # Assign a low value for non-drought days
  
  # Clean up any extreme low values
  problem_rows <- which(data_all$has_drought == 1 & normal_index < -4)
  if (length(problem_rows) > 0) {
    normal_index[problem_rows] <- -4
  }
  
  final_data <- data_all %>% mutate(drought_index = normal_index)
  
  cat("VPD drought index constructed successfully!\n")
  cat("Index range:", round(min(final_data$drought_index), 2), "to", round(max(final_data$drought_index), 2), "\n")
  
  return(list(
    data = final_data,
    vine_model = vine_fit_drought,
    summary_stats = list(
      total_days = nrow(final_data),
      drought_days = sum(final_data$has_drought),
      index_range = range(final_data$drought_index)
    )
  ))
}

# ===============================================================================
# OUTPUT AND VISUALIZATION FUNCTIONS
# ===============================================================================

save_drought_results <- function(final_results, final_config, period_configs) {
  cat("\n===== SAVING VPD DROUGHT RESULTS =====\n")
  
  # Generate output filename based on the full analysis range
  output_filename <- paste0(final_config$output_prefix, "_", final_config$target_start, "_", final_config$target_end, ".csv")
  write.csv(final_results$data, file = output_filename, row.names = FALSE)
  cat("Results saved to:", output_filename, "\n")
  
  # Generate summary filename
  summary_filename <- paste0(final_config$output_prefix, "_summary_", final_config$target_start, "_", final_config$target_end, ".txt")
  
  # Write detailed summary
  sink(summary_filename)
  cat("VPD DROUGHT INDEX CALCULATION SUMMARY\n")
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
  cat("- VPD threshold quantile:", final_config$threshold_quantile, "\n\n")
  
  cat("Results:\n")
  cat("- Total days analyzed:", final_results$summary_stats$total_days, "\n")
  cat("- Days with drought events:", final_results$summary_stats$drought_days, "\n")
  cat("- Index range:", paste(round(final_results$summary_stats$index_range, 2), collapse = " to "), "\n")
  sink()
  
  cat("Summary saved to:", summary_filename, "\n")
}

# ===============================================================================
# MAIN EXECUTION FUNCTION
# ===============================================================================

run_combined_vpd_drought_analysis <- function(master_config, period_configs) {
  cat("===============================================================================\n")
  cat("COMBINED VPD DROUGHT INDEX ANALYSIS\n")
  cat("===============================================================================\n\n")
  
  # Step 1: Load data once
  data <- load_and_prepare_vpd_data(master_config)
  
  # Initialize containers for combined results
  combined_frequency_results <- data.frame()
  combined_severity_results <- data.frame()
  
  # Step 2: Loop through each defined period
  for (i in 1:length(period_configs)) {
    period_name <- names(period_configs)[i]
    period_conf <- period_configs[[i]]
    
    cat(sprintf("\n--- Processing %s (Target: %d-%d, Reference: %d-%d) ---\n\n",
                period_name, period_conf$target_start, period_conf$target_end,
                period_conf$reference_start, period_conf$reference_end))
    
    # Merge master and period-specific configs for the current run
    current_config <- c(master_config, period_conf)
    
    # Validate the configuration for this period
    validate_vpd_config(current_config)
    
    # Calculate daily thresholds based on the specific reference period
    vpd_thresholds <- calculate_daily_vpd_thresholds(data, current_config)
    
    # Calculate drought characteristics for the specific target period
    period_results <- process_all_drought_years(data, vpd_thresholds, current_config)
    
    # Append results to the combined data frames
    combined_frequency_results <- rbind(combined_frequency_results, period_results$frequency)
    combined_severity_results <- rbind(combined_severity_results, period_results$severity)
  }
  
  # Step 3: Combine all period results for final index construction
  combined_all_results <- list(
    frequency = combined_frequency_results,
    severity = combined_severity_results
  )
  
  # Create a final configuration object representing the entire analysis span
  final_config <- master_config
  final_config$target_start <- min(sapply(period_configs, `[[`, "target_start"))
  final_config$target_end <- max(sapply(period_configs, `[[`, "target_end"))
  
  # Step 4: Construct the drought index using the full combined dataset
  final_results <- construct_drought_index(combined_all_results, final_config)
  
  # Step 5: Save the final results and a detailed summary
  save_drought_results(final_results, final_config, period_configs)
  
  cat("\n===============================================================================\n")
  cat("COMBINED VPD DROUGHT ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("===============================================================================\n")
  
  return(final_results)
}

# ===============================================================================
# RUN THE ANALYSIS
# ===============================================================================

# Execute the combined analysis with the configurations defined at the top
drought_results <- run_combined_vpd_drought_analysis(MASTER_CONFIG, PERIOD_CONFIGS)