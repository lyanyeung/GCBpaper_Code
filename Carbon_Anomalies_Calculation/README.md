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
