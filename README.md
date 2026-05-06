# ttt-p1-main-analysis

This repository contains code for the main paper for Phase 1 of Project Track to Treat. The repo is linked to this project on the Open Science Framework (OSF): [https://osf.io/c4e75/](https://osf.io/c4e75/).

# Centralized Data Cleaning

The initial centralized data cleaning code was drafted by [Michael Mullarkey](https://github.com/mcmullarkey).

Lab staff who contributed to Phase I of TTT include former research coordinators Sharon Leong (formerly Chen) and Akash Shroff, and Laura Jans (research assistant at the time).

The data, initial code, and documentation are stored in `jslab/TRACK to TREAT/` on the FSMResFiles server.

## Data

### Raw Qualtrics

Raw baseline and 3-month survey data are stored in `/TRACK to TREAT/Data/Qualtrics Data/Raw Data/`, which contains 18 CSV files obtained from Qualtrics (per Date Modified file metadata, presumably 6 files were obtained 6/16/2020-5/20/2021 and 12 files were obtained on 6/18/2021). The Qualtrics cleaning script appears to focus on the latter 12 files (see below).

### Raw LifePak

Raw EMA data are stored in `/TRACK to TREAT/Data/LifePak Raw Data (Do Not Modify)/`, which contains 10 CSV files obtained from LifePak (per Date Modified file metadata, presumably 8 files were obtained on 4/28/2020 and 2 files were obtained on 9/28/2021). The LifePak cleaning script appears to focus on the 5 files with `NIS` in the filename (see below); the `/TRACK to TREAT/Data/readme_ttt_p1.docx` file also states that files with `NIS` (which it defines as "notification-initiated survey") in the filename are the data to be used.

### Clean

Outputs of the initial data cleaning code are stored in `/TRACK to TREAT/Data/Processed Data/2022 From Michael Mullarkey/`

## Code

TODO: The present repo uses the following scripts from `/TRACK to TREAT/Code/Data Cleaning/old/2022 From Michael Mullarkey/` as a starting point for centralized data cleaning. Given that in general the Qualtrics data seem to have been cleaned before the LifePak data (after which each dataset was deidentified), in this repo the scripts have been numbered in the order to be run.

### `01_ttt_phase1_qualtrics_cleaning.Rmd`

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
- This file, in `/TRACK to TREAT/Code/Data Cleaning/old/2022 From Michael Mullarkey/`, was obtained by Jeremy Eberle from Michael on 10/31/2023. Michael stated that he obtained the file from a Google Drive folder owned by Sharon Chen.

Outputs (though both are commented out) `yb_lsmh_ids_dates.csv` and `cleaned_qualtrics_ttt_phase_1.csv`. Moreover, outputs `cleaned_qualtrics_ttt_phase_1_fixed_220604.csv`, but this does not appear to be used later in data cleaning pipeline (seems later scripts just input `cleaned_qualtrics_ttt_phase_1.csv`)
- Isaac Ahuvia stated that he revised the cleaning script in May 2022 just to keep a variable that had been deleted or something similar, so `cleaned_qualtrics_ttt_phase_1_fixed_220604.csv` may relate to this. Both this script and the separate `ttt_phase1_qualtrics_cleaning_fix.Rmd` (which is not on the present repo but in `/TRACK to TREAT/Code/Data Cleaning/old/2022.06 From Isaac Ahuvia/` and which outputs a CSV file with a different date, `cleaned_qualtrics_ttt_phase_1_fixed_220606.csv`) were last modified on the same date (Date Modified metadata of 6/6/2022).
- Note: `02_ttt_phase1_lifepak_cleaning.Rmd` below inputs `cleaned_qualtrics_ttt_phase_1.csv` and then overwrites it after correcting some participant IDs

### `02_ttt_phase1_lifepak_cleaning.Rmd`

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

### `03_deid_ttt_phase_1.Rmd`

Inputs cleaned LifePak data (`cleaned_lifepak_ttt_phase_1.csv`) and cleaned Qualtrics data (`cleaned_qualtrics_ttt_phase_1.csv`)

Outputs deidentified data (`deid_cleaned_lifepak_ttt_phase_1.csv` and `deid_cleaned_qualtrics_ttt_phase_1.csv`)

## Other Documentation

The following files in `/TRACK to TREAT/` appear relevant to data cleaning

### General

- `/TRACK to TREAT/Data/readme_ttt_p1.docx`
- `/TRACK to TREAT/Data/Processed Data/2022 From Michael Mullarkey/README.rtf`
- `/TRACK to TREAT/Code/Data Cleaning/README_ttt_p1_data_cleaning.docx`
  - Points to [jwe4ec/track-to-treat](https://github.com/jwe4ec/track-to-treat) repo as most recent data cleaning effort

### LifePak

- `/TRACK to TREAT/Data/3TP1_LifePak_Version_IDs.xlsx`
- `/TRACK to TREAT/Data/README info from Laura Jans` folder
  - See contents of this folder for info from Laura Jans re (a) 7 participants who have LifePak data for fewer than the expected number of beeps (see yellow highlights in `2024.04.03 Email with Laura Jans re EMA slider and missing EMA data.pdf`) and (b) whether EMA slider items could be skipped (see orange highlights).
  - For the main TTT paper, these data are treated as missing (see [jwe4ec/ttt-p1-main-analysis](https://github.com/jwe4ec/ttt-p1-main-analysis) repo)

## Issues

- Unable to reproduce clean Qualtrics data
  - As of 12/3/2024, Jeremy can reproduce clean Lifepak Data (`cleaned_lifepak_ttt_phase_1.csv` and `deid_cleaned_lifepak_ttt_phase_1.csv`) per `identical(x, y, FALSE, FALSE, FALSE, FALSE)`
    - Specifically, he can do so using R 4.1.1 (latest version available on 9/28/2021; see below) and the most recent versions of `tidyverse`, `skimr`, `glue`, and `janitor` available on 12/3/2024 (loaded via `library()`). He tried to use the `groundhog` package to load the latest available package versions on 1/7/2022 (date that output files were saved to server; see below) but could not use `groundhog` to load `tidyverse` as `tidyverse` depends on `knitr`, which is "already in use" as it is used to execute Rmd files.
  - However, he cannot reproduce `cleaned_qualtrics_ttt_phase_1.csv`
  - It is unclear what R version and package versions should be used for each script
    - `ttt_phase1_qualtrics_cleaning.Rmd` lists 6/17/2021 as the Date; `ttt_phase1_lifepak_cleaning.Rmd` lists 9/28/2021 as the Date. The output files `cleaned_qualtrics_ttt_phase_1.csv` and `cleaned_lifepak_ttt_phase_1.csv` have Date Modified metadata of 1/7/2022. Thus, the scripts used R and package versions prior to these dates.
    - Note: Michael stated that he cannot guarantee he always used the most up-to-date packages, but he endorsed using these dates as a starting point for determining which R and package versions he used
- Unnecessary packages
  - Although many packages are loaded, only a few appear used by each script (see lists below)
  - Moreover, one loaded package (`doMC`, for parallelization, which does not seem needed) is Unix only and unavailable for Windows
```
# "01_ttt_phase1_qualtrics_cleaning.Rmd" packages: "tidyverse", "glue", "janitor", "fastDummies", "diffdf", "datapasta", "fuzzyjoin"
# "02_ttt_phase1_lifepak_cleaning.Rmd" packages:   "tidyverse", "skimr", "glue", "janitor"
```
- Some packages (`datapasta`, `fuzzyjoin`) are loaded partway through script rather than at top
- One file (`dp5_p1_scoring.csv`) is loaded partway through `ttt_phase1_qualtrics_cleaning.Rmd` rather than at the top
- Hard-coding (e.g., need to use LifePak filenames rather than code below in `ttt_phase1_lifepak_cleaning.Rmd`)
```
files <- list.files(pattern = "*.csv")
lifepak_files <- files[c(1:5)]
```


## Network Analyses and Prediction Models

**TODO: Jeremy to update this section**





### Data

Analyses input cleaned and deidentified LifePak data in `deid_cleaned_lifepak_ttt_phase_1.csv`

### Code

Initial network analyses were run by [Sebastian Castro-Alvarez](https://github.com/secastroal) and [Laura
Bringmann](https://github.com/LauraBringmann) using code in TODO. Analyses were revised by Josip Razum and [Jeremy Eberle](https://github.com/jwe4ec).

**TODO: Jeremy to document scripts**





## Other Machine Learning

**TODO: Jeremy to update this section**





### Data

Analyses input cleaned (but not deidentified) LifePak data in `cleaned_lifepak_ttt_phase_1.csv` and cleaned (but not deidentified) Qualtrics data in `cleaned_qualtrics_ttt_phase_1.csv`.
- TODO: Consider revising code to import deidentified data instead)

### Code

Initial machine learning analyses were run by [Michael Mullarkey](https://github.com/mcmullarkey) using
code located at TODO. Analyses were revised by [Yama Chang](https://github.com/yamachang).

#### Michael's

**TODO: Compare code on server with code on Michael's [ttt-main-analyses](https://github.com/mcmullarkey/ttt_main_analyses) repo to help determine if code on server is Michael's original code**
- See `jslab\TRACK to TREAT\Code\Primary Analyses`

#### Yama's

**TODO: Process and describe Yama's latest code in folder below**
- See `jslab\TRACK to TREAT\Code\Primary Analyses\2024.08.01 From Yama Chang`