## Extreme event index calculation

The folder `Extreme_Event_Index_Calculation/` contains R scripts for calculating heatwave and atmospheric drought indices for FLUXNET sites.

## Files

- `Extreme_Event_Index_Calculation/Heatwave_single_site_example.R`: single-site example script for calculating the heatwave index.
- `Extreme_Event_Index_Calculation/Heatwave_batch_run.R`: batch-processing script for calculating the heatwave index across multiple sites.
- `Extreme_Event_Index_Calculation/AtmosphericDrought_single_site_example.R`: single-site example script for calculating the atmospheric drought index.
- `Extreme_Event_Index_Calculation/AtmosphericDrought_batch_run.R`: batch-processing script for calculating the atmospheric drought index across multiple sites.
- `Extreme_Event_Index_Calculation/folder_list.csv`: list of site folders used by the batch-processing scripts.

## Input data

The scripts require site-level ERA5/ERA5-Land meteorological input CSV files.

The raw ERA5/ERA5-Land input data are not included in this repository because of file size.

Required columns for heatwave index calculation:

- `DATETIME`
- `temperature_2m_C`

Required columns for atmospheric drought index calculation:

- `DATETIME`
- `VPD_kPa`

The input files are assumed to follow this naming pattern:

`SITE_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv`

Example:

`DK-Sor_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv`

## Folder list

The batch scripts use `folder_list.csv`, which should contain a column named `folder`.

Example:

| folder |
|---|
| DK_Sor_1996_2020 |
| FR_Pue_2000_2020 |
| IT_Ro1_2002_2020 |

## Expected local data structure

For heatwave batch processing, organise the local input data as:

- `~/data/heatwave_inputs/folder_list.csv`
- `~/data/heatwave_inputs/DK_Sor_1996_2020/DK-Sor_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv`
- `~/data/heatwave_inputs/FR_Pue_2000_2020/FR-Pue_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv`

For atmospheric drought batch processing, organise the local input data as:

- `~/data/atmospheric_drought_inputs/folder_list.csv`
- `~/data/atmospheric_drought_inputs/DK_Sor_1996_2020/DK-Sor_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv`
- `~/data/atmospheric_drought_inputs/FR_Pue_2000_2020/FR-Pue_ERA5_Land_stitched_1960-12-16_to_2025-10-10.csv`

## Required R packages

Install the required R packages before running the scripts:

`install.packages(c("tidyverse", "lubridate", "ggplot2", "gridExtra", "dplyr", "corrplot", "FactoMineR", "factoextra", "viridis", "patchwork", "moments", "RColorBrewer", "scales", "MASS", "scatterplot3d", "copula", "VineCopula", "kdecopula", "purrr"))`

## How to run

### Single-site heatwave calculation

Open:

`Extreme_Event_Index_Calculation/Heatwave_single_site_example.R`

Edit the `MASTER_CONFIG` section at the beginning of the script, especially:

- `working_dir`
- `input_file`
- `output_prefix`

Then run in R:

`source("Extreme_Event_Index_Calculation/Heatwave_single_site_example.R")`

### Batch heatwave calculation

Open:

`Extreme_Event_Index_Calculation/Heatwave_batch_run.R`

Edit the `BATCH_CONFIG` section at the beginning of the script, especially:

- `base_dir`
- `folder_list_file`
- `output_base_dir`
- `source_file`

Then run in R:

`source("Extreme_Event_Index_Calculation/Heatwave_batch_run.R")`

### Single-site atmospheric drought calculation

Open:

`Extreme_Event_Index_Calculation/AtmosphericDrought_single_site_example.R`

Edit the `MASTER_CONFIG` section at the beginning of the script, especially:

- `working_dir`
- `input_file`
- `output_prefix`

Then run in R:

`source("Extreme_Event_Index_Calculation/AtmosphericDrought_single_site_example.R")`

### Batch atmospheric drought calculation

Open:

`Extreme_Event_Index_Calculation/AtmosphericDrought_batch_run.R`

Edit the `BATCH_CONFIG` section at the beginning of the script, especially:

- `base_dir`
- `folder_list_file`
- `output_base_dir`
- `source_file`

Then run in R:

`source("Extreme_Event_Index_Calculation/AtmosphericDrought_batch_run.R")`

## Output

The scripts generate CSV output files containing the calculated heatwave or atmospheric drought indices.

For batch processing, outputs are saved to the directory specified by `output_base_dir` in the corresponding `BATCH_CONFIG`.

## Notes

Users should modify the file paths in `MASTER_CONFIG` or `BATCH_CONFIG` before running the scripts.

The raw ERA5/ERA5-Land data are not included in this repository because of file size.

The batch scripts assume that each site folder contains an ERA5/ERA5-Land input CSV file using the site code naming pattern described above.
