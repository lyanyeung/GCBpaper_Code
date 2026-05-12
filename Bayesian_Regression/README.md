# Bayesian_Regression

This repository contains the R scripts and input data used for Bayesian meta-regression analyses of carbon anomaly responses under heatwave and VPD drought scenarios.

The analyses use Gaussian Bayesian meta-regression models to evaluate how site-level environmental and ecological predictors modulate carbon anomalies, represented as Z-scores.

## Repository structure

```text
Bayesian_Regression/
├── README.md
├── run_Bayesian_MetaRegression_H_Gaussian.R
├── run_Bayesian_MetaRegression_VPD_Gaussian.R
├── data/
│   ├── Folder_List_All_North.csv
│   ├── All_Sites_ZScore_H_Statistics_Summary.csv
│   └── All_Sites_ZScore_VPD_Statistics_Summary.csv
└── outputs/
    └── .gitkeep
```

## Scripts

This repository contains two main R scripts.

### 1. Heatwave analysis

```text
run_Bayesian_MetaRegression_H_Gaussian.R
```

This script runs the Gaussian Bayesian meta-regression analysis for heatwave-related scenarios:

1. `HW_ge_2`: heatwave condition, defined as H >= 2
2. `HW_ge_2_VPD_lt_2`: pure heatwave condition, defined as H >= 2 and VPD < 2
3. `HW_ge_2_VPD_ge_2`: compound heatwave and VPD drought condition, defined as H >= 2 and VPD >= 2

The script uses the following input files:

```text
data/All_Sites_ZScore_H_Statistics_Summary.csv
data/Folder_List_All_North.csv
```

Outputs are saved to:

```text
outputs/Gaussian_Arbitrary_Combinations_North_With_LAT_v2/
```

### 2. VPD drought analysis

```text
run_Bayesian_MetaRegression_VPD_Gaussian.R
```

This script runs the Gaussian Bayesian meta-regression analysis for VPD drought-related scenarios:

1. `VPD_ge_2`: VPD drought condition, defined as VPD >= 2
2. `HW_lt_2_VPD_ge_2`: pure VPD drought condition, defined as VPD >= 2 and H < 2
3. `HW_ge_2_VPD_ge_2`: compound VPD drought and heatwave condition, defined as VPD >= 2 and H >= 2

The script uses the following input files:

```text
data/All_Sites_ZScore_VPD_Statistics_Summary.csv
data/Folder_List_All_North.csv
```

Outputs are saved to:

```text
outputs/Gaussian_Arbitrary_Combinations_VPD_North_With_LAT/
```

## Input data

The `data/` folder contains three required input files.

### 1. Site-level metadata

```text
Folder_List_All_North.csv
```

This file contains site-level environmental and ecological predictors used in the Bayesian meta-regression models.

Required columns include:

```text
folder
IGBP
canopy_height
Mean_Annual_Precip
Growing_LAI
Fpar
MAT
Weighted_AWM
LATITUDE
```

The scripts use the following predictors:

| Variable in model | Source column |
|---|---|
| CH | canopy_height |
| MAP | Mean_Annual_Precip |
| LAI | Growing_LAI |
| Fpar | Fpar |
| MAT | MAT |
| AWM | Weighted_AWM |
| LAT | LATITUDE |

### 2. Heatwave Z-score summary file

```text
All_Sites_ZScore_H_Statistics_Summary.csv
```

This file contains site-level Z-score response summaries for heatwave-related scenarios.

Required columns include:

```text
site_code
scenario
n
mean_zscore
zscore_se
p_value
```

The Bayesian models use:

| Column | Description |
|---|---|
| site_code | Site identifier |
| scenario | Scenario name |
| n | Sample size |
| mean_zscore | Response variable |
| zscore_se | Standard error of the response variable |
| p_value | Site-level significance value |

### 3. VPD Z-score summary file

```text
All_Sites_ZScore_VPD_Statistics_Summary.csv
```

This file contains site-level Z-score response summaries for VPD drought-related scenarios.

Required columns include:

```text
site_code
scenario
n
mean_zscore
zscore_se
p_value
```

The Bayesian models use:

| Column | Description |
|---|---|
| site_code | Site identifier |
| scenario | Scenario name |
| n | Sample size |
| mean_zscore | Response variable |
| zscore_se | Standard error of the response variable |
| p_value | Site-level significance value |

## Model description

For each scenario, the scripts fit Gaussian Bayesian meta-regression models using the `brms` package.

The response variable is:

```text
mean_zscore
```

with uncertainty represented by:

```text
zscore_se
```

In the scripts, these variables are renamed as:

```r
beta_hw = mean_zscore
se_beta = zscore_se
```

The fitted model form is:

```r
beta_hw | se(se_beta, sigma = TRUE) ~ predictors
```

The candidate predictors are:

```text
CH_std
MAP_std
LAI_std
Fpar_std
MAT_std
AWM_std
LAT_std
```

All predictors are standardized before modelling.

The scripts automatically generate all predictor combinations from one-variable models up to four-variable models. With seven predictors, this results in all possible combinations of one, two, three, and four predictors for each scenario.

## Requirements

The analyses were written in R and require the following R packages:

```r
tidyverse
ggplot2
brms
patchwork
posterior
ggpmisc
utils
```

Before running the scripts, please install the required packages if they are not already installed:

```r
install.packages(c(
  "tidyverse",
  "ggplot2",
  "brms",
  "patchwork",
  "posterior",
  "ggpmisc"
))
```

The `utils` package is included with base R.

Because the models are fitted using `brms`, a working Stan backend is required. Depending on the local R setup, users may need to configure either `rstan` or `cmdstanr`.

## How to run

Clone or download this repository, then open R or RStudio.

Set the working directory to the root folder of the repository:

```r
setwd("path/to/Bayesian_Regression")
```

For example, the working directory should be the folder that contains:

```text
run_Bayesian_MetaRegression_H_Gaussian.R
run_Bayesian_MetaRegression_VPD_Gaussian.R
data/
outputs/
```

### Run the heatwave analysis

```r
source("run_Bayesian_MetaRegression_H_Gaussian.R")
```

### Run the VPD drought analysis

```r
source("run_Bayesian_MetaRegression_VPD_Gaussian.R")
```

The scripts use relative paths. Therefore, no user-specific absolute paths are required.

## Output files

Each script creates an output folder inside `outputs/`.

The heatwave script saves results to:

```text
outputs/Gaussian_Arbitrary_Combinations_North_With_LAT_v2/
```

The VPD drought script saves results to:

```text
outputs/Gaussian_Arbitrary_Combinations_VPD_North_With_LAT/
```

Each analysis produces:

1. A detailed model summary table as a `.csv` file
2. Scatter plots for each predictor as `.png` files
3. A forest plot for univariate models as a `.png` file
4. A saved R workspace as a `.RData` file

Example output files from the heatwave analysis include:

```text
AllCombos_H_ZScore_Gaussian_WithLAT_results_summary.csv
AllCombos_H_ZScore_Gaussian_WithLAT_Forest_Univariate.png
AllCombos_H_ZScore_Gaussian_WithLAT_workspace.RData
```

Example output files from the VPD drought analysis include:

```text
AllCombos_VPD_ZScore_Gaussian_WithLAT_results_summary.csv
AllCombos_VPD_ZScore_Gaussian_WithLAT_Forest_Univariate.png
AllCombos_VPD_ZScore_Gaussian_WithLAT_workspace.RData
```

The saved workspace contains:

```text
results_summary
all_scenario_data
all_models_objects
```

The object `all_models_objects` contains the fitted `brms` model objects.

## Notes for reviewers

The scripts are designed to be run directly from the root directory of this repository.

All input files should remain in the `data/` folder. The scripts use relative paths and do not require local absolute paths such as `D:/...` or `E:/...`.

The full Bayesian model fitting process may take a long time because each script fits all predictor combinations up to four variables for each scenario. Runtime will depend on the machine, number of CPU cores, and Stan configuration.

By default, each model is fitted with:

```r
iter = 4000
warmup = 1000
chains = 4
cores = 4
seed = 123
```

For reproducibility, the random seed is fixed in the model fitting step.

## Recommended `.gitignore`

Large output files and fitted model workspaces can be excluded from GitHub if desired. A recommended `.gitignore` is:

```gitignore
outputs/*
!outputs/.gitkeep

*.RData
*.RDS
.Rhistory
.RData
.Rproj.user/
```

If the fitted model objects are required for review, the `.RData` files can be shared separately because they may be large.

## License

Please refer to the project or manuscript documentation for data usage and licensing information.
