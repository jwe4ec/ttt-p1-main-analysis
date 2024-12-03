# ttt-p1-main-analysis

This repository contains analysis code for this project on the Open Science Framework (OSF): [https://osf.io/c4e75/](https://osf.io/c4e75/).

## Table of Contents

- [Project Overview](#project-overview)
- [Centralized Data Cleaning](#centralized-data-cleaning)
  - [Data](#data)
  - [Code](#code)
  - [Other Documentation](#other-documentation)
  - [TODOs](#todos)
- [Network Analyses](#network-analyses)
- [Machine Learning](#machine-learning)

## Project Overview

Phase I of Project TRACK to TREAT aims to use parameters from network models estimated from ecological momentary assessment (EMA) data to predict 3-month changes in depression symptoms and related constructs in depressed adolescents. Phase I, an observational study, consisted of a baseline Qualtrics survey, 21 days of EMA (5 pings per day), and another Qualtrics survey 3 months later.

## Centralized Data Cleaning

Data, initial code, and documentation relevant to centralized data cleaning for Phase I of Project TRACK to TREAT (TTT) are stored in the `MSS/Schleider_Lab/jslab/TRACK to TREAT` folder on the [FSMResFiles](https://www.feinberg.northwestern.edu/it/services/server-storage-and-data/research-data-storage.html) server at [Northwestern University Feinberg School of Medicine](https://www.feinberg.northwestern.edu/).

The initial centralized data cleaning code was drafted by [Michael Mullarkey](https://github.com/mcmullarkey). The Centralized Data Cleaning section of the present repo houses [Jeremy Eberle](https://github.com/jwe4ec) and [Isaac Ahuvia](https://github.com/isaacahuvia)'s attempts to improve the code for greater reproducibility. For centralized data cleaning for Phase II of TTT, see the separate repo [ttt-p2-cleaning](https://github.com/jwe4ec/ttt-p2-cleaning).

Former lab staff who contributed to Phase I of TTT include Sharon Chen (research coordinator at the time) and Laura Jans (research assistant at the time).

### Data

#### Raw

##### From Qualtrics

Raw baseline and 3-month survey data are stored in the `/TRACK to TREAT/Data/Qualtrics Data/Raw Data` folder, which contains 18 CSV files obtained from Qualtrics (per Date Modified file metadata, presumably 6 files were obtained 6/16/20-5/20/21 and 12 files were obtained on 6/18/21). The Qualtrics cleaning script appears to focus on the latter 12 files (see below).

##### From LifePak

Raw EMA data are stored in the `/TRACK to TREAT/Data/LifePak Raw Data (Do Not Modify)` folder, which contains 10 CSV files obtained from LifePak (per Date Modified file metadata, presumably 8 files were obtained on 4/28/20 and 2 files were obtained on 9/28/21). The LifePak cleaning script appears to focus on the 5 files with `NIS` in the filename (see below); the `/TRACK to TREAT/Data/readme_ttt_p1.docx` file also states that files with `NIS` (which it defines as "notification-initiated survey") in the filename are the data to be used.

#### Clean

Outputs of the initial data cleaning code are stored in the `/TRACK to TREAT/Data/Processed Data` folder

### Code

The present repo uses the following scripts from the `/TRACK to TREAT/Code/Cleaning Data` folder as a starting point for centralized data cleaning. Given that in general the Qualtrics data seem to have been cleaned before the LifePak data (after which each dataset was deidentified), in this repo the scripts have been numbered in the order to be run.

#### `01_ttt_phase1_qualtrics_cleaning.Rmd`

Inputs the following 12 raw CSV files (out of the 18 from Qualtrics)
```
# "dp5_b_child_p1_numeric.csv"
# "dp5_b_child_p1_choice_text.csv"
# "dp5_b_child_remote_p1_numeric.csv"
# "dp5_b_child_remote_p1_choice_text.csv"
# "dp5_b_parent_p1_numeric.csv"
# "dp5_b_parent_p1_choice_text.csv"
# "dp5_b_parent_remote_p1_numeric.csv"
# "dp5_b_parent_remote_p1_choice_text.csv"
# "dp5_3m_child_p1_numeric.csv"
# "dp5_3m_child_p1_choice_text.csv"
# "dp5_3m_parent_p1_numeric.csv"
# "dp5_3m_parent_p1_choice_text.csv"
```

Also inputs `dp5_p1_scoring.csv`
- This file, in `/TRACK to TREAT/Code/Cleaning Data`, was obtained by Jeremy Eberle from Michael Mullarkey on 10/31/23. Michael stated that he obtained the file from a Google Drive folder owned by Sharon Chen.

Outputs (though both are commented out) `yb_lsmh_ids_dates.csv` and `cleaned_qualtrics_ttt_phase_1.csv`. Moreover, outputs `cleaned_qualtrics_ttt_phase_1_fixed_220604.csv`, but this does not appear to be used later in data cleaning pipeline (seems later scripts just input `cleaned_qualtrics_ttt_phase_1.csv`)
- Isaac Ahuvia stated that he revised the cleaning script in May 2022 just to keep a variable that had been deleted or something similar, so `cleaned_qualtrics_ttt_phase_1_fixed_220604.csv` may relate to this. Both this script and a separate script `ttt_phase1_qualtrics_cleaning_fix.Rmd` (which is not on the present repo and which outputs a CSV file with a different date, `cleaned_qualtrics_ttt_phase_1_fixed_220606.csv`) were last modified on the same date (Date Modified metadata of 6/6/22).
- Note: `02_ttt_phase1_lifepak_cleaning.Rmd` below inputs `cleaned_qualtrics_ttt_phase_1.csv` and then overwrites it after correcting some participant IDs

#### `02_ttt_phase1_lifepak_cleaning.Rmd`

Inputs the following 5 raw CSV files (out of the 10 from LifePak)
```
# "3T_P1_V1_NIS_2020_Mar_02.csv"
# "3T_P1_V2_NIS_2020_Mar_13.csv"
# "3T_P1_V2_NIS_21200_958251_Download2.csv"
# "3T_P1_V2_NIS_21200_958251_Download3.csv"
# "3T_P1_V4_NIS.csv"
```

Also inputs `cleaned_qualtrics_ttt_phase_1.csv` (presumably originally from `01_ttt_phase1_qualtrics_cleaning.Rmd`)

Outputs `cleaned_combined_qualtrics_lifepak_ttt_phase_1.csv` (does not appear to be used later in data cleaning pipeline), `cleaned_lifepak_ttt_phase_1.csv`, and `cleaned_qualtrics_ttt_phase_1.csv`

#### `03_deid_ttt_phase_1.Rmd`

Inputs cleaned LifePak data (`cleaned_lifepak_ttt_phase_1.csv`) and cleaned Qualtrics data (`cleaned_qualtrics_ttt_phase_1.csv`)

Outputs deidentified data (`deid_cleaned_lifepak_ttt_phase_1.csv` and `deid_cleaned_qualtrics_ttt_phase_1.csv`)

### Other Documentation

The following files in the `MSS/Schleider_Lab/jslab/TRACK to TREAT` folder appear relevant to data cleaning

#### General

- `TRACK to TREAT/Data/readme_ttt_p1.docx`
- `TRACK to TREAT/Data/Processed Data/README.rtf`

#### LifePak

- `TRACK to TREAT/Data/3TP1_LifePak_Version_IDs.xlsx`
- `TRACK to TREAT/Data/README info from Laura Jans` folder
  - See contents of this folder for info from Laura Jans re (a) 7 participants who have LifePak data for fewer than the expected number of beeps (see yellow highlights in `2024.04.03 Email with Laura Jans re EMA slider and missing EMA data.pdf`) and (b) whether EMA slider items could be skipped (see orange highlights).
  - For the main TTT paper, these data are treated as missing (see `02_networks/code/02_further_clean_data_align_obs.R` in the present repo)

### TODOs

#### General

- TODO: As of 12/3/24, Jeremy can reproduce `cleaned_lifepak_ttt_phase_1.csv` (and `deid_cleaned_lifepak_ttt_phase_1.csv`) per `identical(x, y, FALSE, FALSE, FALSE, FALSE)`. However, he cannot reproduce `cleaned_qualtrics_ttt_phase_1.csv`.
  - Specifically, he can reproduce the clean LifePak data when using R 4.1.1 (latest version available on 9/28/21; see below) and the most recent versions of `tidyverse`, `skimr`, `glue`, and `janitor` available on 12/3/24 (loaded via `library()`). He tried to use the `groundhog` package to load the latest available package versions on 1/7/22 (date that output files were saved to server; see below) but could not use `groundhog` to load `tidyverse` as `tidyverse` depends on `knitr`, which is "already in use" as it is used to execute Rmd files.

#### Specific

- TODO: Determine what R version and package versions should be used for each script
  - `01_ttt_phase1_qualtrics_cleaning.Rmd` lists 6/17/2021 as the Date; `02_ttt_phase1_lifepak_cleaning.Rmd` lists 9/28/2021 as the Date. The output files `cleaned_qualtrics_ttt_phase_1.csv` and `cleaned_lifepak_ttt_phase_1.csv` have Date Modified metadata of 1/7/22. Thus, the scripts used R and package versions prior to these dates.
  - Note: Michael stated that he cannot guarantee he always used the most up-to-date packages, but he endorsed using these dates as a starting point for determining which R and package versions he used
- TODO: Determine what packages are needed, load only those, and load all needed packages at top of script
  - Although many packages are loaded, only a few appear used by each script (see lists below)
  - Moreover, one of the loaded packages (`doMC`, for parallelization) is Unix only and unavailable for Windows
  - Some packages (`datapasta`, `fuzzyjoin`) are loaded partway through script rather than at top
  - Find a way to control the version of `knitr` (needed for Rmd files) or do not use Rmd files (see issue above)
```
# "01_ttt_phase1_qualtrics_cleaning.Rmd" packages: "tidyverse", "glue", "janitor", "fastDummies", "diffdf", "datapasta", "fuzzyjoin"
# "02_ttt_phase1_lifepak_cleaning.Rmd" packages:   "tidyverse", "skimr", "glue", "janitor"
```
- TODO: Remove parallelization, as it does not seem needed
- TODO: Treating `.` as the parent folder for the present repo, create local raw data folders (`./data/raw/qualtrics` and `./data/raw/lifepak`), put raw CSV files in those folders, and use relative file paths to load raw data and output clean data (vs. storing raw data, clean data, and code in same folder)--see below for example. We can then describe the directory structure in this README.
```
example_raw_table <- read.csv("./data/raw/qualtrics/example_raw_table.csv")

clean_path <- "./data/clean/"
dir.create(clean_path)
write.csv(example_clean_table, paste0(clean_path, "example_clean_table.csv"))
```
- TODO: Load `dp5_p1_scoring.csv` at top of `01_ttt_phase1_qualtrics_cleaning.Rmd` (vs. partway through script)
- TODO: Clearly reflect what `01_ttt_phase1_qualtrics_cleaning.Rmd` should output (see description of its outputs above for various issues)
- TODO: Remove extraneous code/comments
- TODO: Avoid hard-coding practices
- TODO: Put deidentified clean data in `./data/clean` folder on [OSF project](https://osf.io/c4e75/) linked to the present repo





## Network Analyses

### Data

Analyses input cleaned and deidentified LifePak data in `deid_cleaned_lifepak_ttt_phase_1.csv`

### Code

Initial network analyses were run by [Sebastian Castro-Alvarez](https://github.com/secastroal) and [Laura
Bringmann](https://github.com/LauraBringmann) using code in the `from_sebastian/` folder. Analyses were revised by Josip Razum and [Jeremy Eberle](https://github.com/jwe4ec).

**TODO: Jeremy to document scripts**





## Machine Learning

**TODO: Jeremy to update this section**





### Data

Analyses input cleaned (but not deidentified) LifePak data in `cleaned_lifepak_ttt_phase_1.csv` and cleaned (but not deidentified) Qualtrics data in `cleaned_qualtrics_ttt_phase_1.csv`.
- TODO: Consider updating `ttt_phase_1_main_analyses_08232023_final.Rmd` to input deidentified data in `deid_cleaned_lifepak_ttt_phase_1.csv` and `deid_cleaned_qualtrics_ttt_phase_1.csv`

### Code

Initial machine learning analyses were run by [Michael Mullarkey](https://github.com/mcmullarkey) using
code located at TODO. Analyses were revised by [Yama Chang](https://github.com/yamachang).

#### `ttt_phase_1_main_analyses_08232023_final.Rmd` (in progress and not yet on repo)

Inputs `cleaned_lifepak_ttt_phase_1.csv`, `cleaned_qualtrics_ttt_phase_1.csv`, and `extracted_features_r.csv`
- TODO: Where is `extracted_features_r.csv` created?

Outputs (though all are commented out) `train_long_ema.csv`, `cdi_data_init.csv`, and `analytical_base_table.csv`

#### `Extract feature.ipynb` (in progress and not yet on repo)

Inputs `train_long_ema.csv`, `test_long_ema.csv`, and `extracted_features_r.csv`
- TODO: Where are `test_long_ema.csv` and `extracted_features_r.csv` created?

Outputs `extracted_features_train_0908.csv` and `extracted_features_test_0908.csv`