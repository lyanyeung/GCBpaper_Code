# Carbon Anomalies Calculation

This repository contains R scripts for calculating carbon flux anomalies and carbon uptake intervals from site-level flux data.

## Files

- `CarbonAnomaly_Batch.R`  
  Calculates daily anomalies for carbon flux variables, including NEE, GPP, and RECO.

- `Carbonuptakeinterval_Batch.R`  
  Calculates site-specific climatological seasonal cycles and identifies carbon uptake intervals.

- `Folder_List_All.csv`  
  Provides the list of site folders and site-level metadata.

## Input data structure

The scripts expect the data to be organised as:

```text
data/
├── Folder_List_All.csv
├── AT_Neu_2002_2020/
│   └── AT_Neu.csv
├── AU_How_2002_2024/
│   └── AU_How.csv
└── ...
```

## Required input columns

Each site CSV file should include:

- `TIMESTAMP`
- `NEE_VUT_REF`
- `GPP_NT_VUT_REF`
- `RECO_NT_VUT_REF`

For carbon anomaly calculation, the following quality-control column is also required:

- `NEE_VUT_REF_QC`

Missing values coded as `-9999` are treated as `NA`.

## Required R packages

The scripts require the following R packages:

```r
install.packages(c("tidyverse", "lubridate", "ggplot2", "mgcv", "rlang", "qpdf"))
```

## How to run
update the input and output paths in the R scripts, then run:
```r
source("CarbonAnomaly_Batch.R")
source("Carbonuptakeinterval_Batch.R")
