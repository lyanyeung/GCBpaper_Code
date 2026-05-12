Carbon Anomalies Calculation

This folder contains the R scripts used to calculate daily carbon flux anomalies and carbon uptake intervals for multiple eddy-covariance sites.

Files
CarbonAnomaly_Batch.R: calculates daily anomalies for NEE, GPP, and ecosystem respiration.
Carbonuptakeinterval_Batch.R: identifies carbon uptake intervals based on climatological NEE < 0.
Folder_List_All.csv: lists all site folders used in the batch processing.
Input Data

Raw site-level flux data are not included in this repository due to data sharing restrictions.

Users should prepare the input data using the following structure:

base_dir/
Folder_List_All.csv
AU_How_2001_2014/
AU_How.csv
BE_Bra_1997_2014/
BE_Bra.csv

Each site folder should be named as:

SITE_YYYY_YYYY

For example:

AU_How_2001_2014

The years in the folder name are used as the baseline, analysis, detrending, and export periods.

Each site folder should contain one CSV file named using the site code, for example:

AU_How_2001_2014/AU_How.csv

Required Columns

Each site-level CSV file should contain:

TIMESTAMP
NEE_VUT_REF
NEE_VUT_REF_QC
GPP_NT_VUT_REF
RECO_NT_VUT_REF

TIMESTAMP should be in YYYYMMDD format.

Missing values coded as -9999 are treated as NA.

Carbon Anomaly Calculation

Run:

source("CarbonAnomaly_Batch.R")

This script calculates daily anomalies for:

NEE_VUT_REF
GPP_NT_VUT_REF
RECO_NT_VUT_REF

The main steps are:

Read each site file listed in Folder_List_All.csv.
Apply quality-control filtering using NEE_VUT_REF_QC >= 0.7.
Calculate climatological seasonal cycles using a month-day based 31-day moving window.
Calculate residuals as observed value minus climatological mean.
Fit and remove a linear temporal trend.
Calculate daily anomalies and standardized z-scores.

For each site, the main output is:

carbon_anomalies_QC_filtered_SITE.csv

The script also generates PDF plots showing residual trends and anomaly time series.

Carbon Uptake Interval Calculation

Run:

source("Carbonuptakeinterval_Batch.R")

This script identifies carbon uptake intervals based on climatological NEE.

Carbon uptake intervals are defined as contiguous periods when climatological NEE < 0.

Negative NEE indicates net ecosystem carbon uptake.

The main output is:

All_Sites_Carbon_Uptake_Intervals.csv

This output includes:

site_code
start_doy
end_doy
start_monthday
end_monthday

The script also creates climatology plots for each site and combines them into:

All_Sites_Combined_Climatology.pdf

Month-day Climatology Method

Both scripts use a month-day based climatology workflow.

Each calendar day is represented by a month-day key, such as 01-01, 01-02, and 12-31.

For each target calendar day, climatology is calculated using a centred 31-day moving window across all baseline years.

February 29 is included only when leap-year data are available.

Configuration

Before running the scripts, users should edit the path settings near the top of each R script.

For example:

base_dir <- "path/to/your/data_directory"

output_base_dir <- "path/to/your/output_directory"

base_dir should point to the folder containing Folder_List_All.csv and all site folders.

output_base_dir is the folder where results will be saved.

Folder List File

Folder_List_All.csv should contain at least one column named folder.

Example:

folder
AU_How_2001_2014
BE_Bra_1997_2014
CH_Dav_1997_2014

Running Order

The two scripts can be run independently because each script recalculates its own climatological baseline from the raw site-level data.

Recommended order:

CarbonAnomaly_Batch.R
Carbonuptakeinterval_Batch.R
R Packages

The scripts require the following R packages:

tidyverse
lubridate
ggplot2
mgcv
rlang
qpdf

qpdf is only required by Carbonuptakeinterval_Batch.R.

Data Availability

Raw flux data are not included in this repository due to data sharing restrictions.

Users should obtain the required site-level flux data separately and place them in the expected folder structure before running the scripts.
