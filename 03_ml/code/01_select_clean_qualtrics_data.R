# ---------------------------------------------------------------------------- #
# Select Clean Qualtrics Data for Prediction Models -----
# Author: Jeremy W. Eberle
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# ---------------------------------------------------------------------------- #
# Store working directory, check correct R version, load packages ----
# ---------------------------------------------------------------------------- #

# Store working directory

wd_dir <- getwd()

# Load custom functions

source("./02_networks/code/01_define_functions.R")

# Check correct R version, load groundhog package, and specify groundhog_day

groundhog_day <- version_control()

# No packages loaded

# ---------------------------------------------------------------------------- #
# Import clean Qualtrics data  ----
# ---------------------------------------------------------------------------- #

# TODO: Finalize file organization after Phase I data cleaning is complete





# Note: Must be connected to VPN to access ResFiles directory below

clean_qualtrics_dat_path <- "R:/MSS/Schleider_Lab/jslab/TRACK to TREAT/Data/Clean Data (Isaac)/"

y_qualtrics_dat <- readRDS(paste0(clean_qualtrics_dat_path, "Phase 1 Youth Qualtrics Data.RDS"))
p_qualtrics_dat <- readRDS(paste0(clean_qualtrics_dat_path, "Phase 1 Parent Qualtrics Data.RDS"))

# ---------------------------------------------------------------------------- #
# Define columns of interest ----
# ---------------------------------------------------------------------------- #

# Define metadata columns

y_meta_cols <- c("lsmh_id", "lifepak_id", "yb_administration", "yb_date", "y3m_date")

p_meta_cols <- c("lsmh_id", "pb_administration", "pb_date", "p3m_date")

all(y_meta_cols %in% names(y_qualtrics_dat))
all(p_meta_cols %in% names(p_qualtrics_dat))

# TODO (Check factor structure): Define relevant outcome columns at baseline and 3 months. 
# Note: Reverse-coded columns have already been unreversed (as part of data cleaning).





# TODO (get all items/scores for a given measure and then exclude irrelevant ones): Youth-rated





y_cdi_cols     <- c(paste0("yb_cdi_",  1:28),
                    paste0("y3m_cdi_", 1:28),
                    "yb_cdi_mean", "y3m_cdi_mean")

y_bhs_cols     <- c(paste0("yb_bhs_",  1:4),
                    paste0("y3m_bhs_", 1:4),
                    "yb_bhs_mean", "y3m_bhs_mean")

y_pcsc_cols    <- c(paste0("yb_pcsc_",  1:24),
                    paste0("y3m_pcsc_", 1:24),
                    "yb_pcsc_mean", "y3m_pcsc_mean")

y_bads_ac_cols <- c(paste0("yb_bads_",  c(3, 4, 5, 7, 11, 12, 23)),
                    paste0("y3m_bads_", c(3, 4, 5, 7, 11, 12, 23)),
                    "yb_bads_ac_mean", "y3m_bads_ac_mean")

y_bads_ar_cols <- c(paste0("yb_bads_",  c(8, 9, 10, 13, 14, 15, 24, 25)),
                    paste0("y3m_bads_", c(8, 9, 10, 13, 14, 15, 24, 25)),
                    "yb_bads_ar_mean", "y3m_bads_ar_mean")

y_outcome_cols <- c(y_cdi_cols, y_bhs_cols, y_pcsc_cols, y_bads_ac_cols, y_bads_ar_cols)

all(y_outcome_cols %in% names(y_qualtrics_dat))

# TODO (get all demographics items and then exclude irrelevant ones): Parent-rated





p_cdi_cols     <- c(paste0("pb_cdi_",  1:17),
                    paste0("p3m_cdi_", 1:17),
                    "pb_cdi_mean", "p3m_cdi_mean")

# TODO: Add parent demographics once cleaned (at minimum, parent age, sex, gender, 
# race [if assessed], ethnicity, single- vs. co-parent, depression symptom severity).
# Also, further cleaning of "pb_childgender" and "pb_school" are in progress.





p_dem_cols     <- c("pb_childage", "pb_childsex", "pb_childgender", "pb_childethnicity", 
                    "pb_grade", "pb_school", "pb_income", "pb_dependent")

all(c(p_cdi_cols, p_dem_cols) %in% names(p_qualtrics_dat))

# ---------------------------------------------------------------------------- #
# Select columns of interest ----
# ---------------------------------------------------------------------------- #

y_qualtrics_dat_sel <- y_qualtrics_dat[, c(y_meta_cols, y_outcome_cols)]
p_qualtrics_dat_sel <- p_qualtrics_dat[, c(p_meta_cols, p_cdi_cols, p_dem_cols)]

# ---------------------------------------------------------------------------- #
# Save data ----
# ---------------------------------------------------------------------------- #

selected_clean_dat_path <- "./03_ml/data/selected_clean/"

dir.create(selected_clean_dat_path, recursive = TRUE)

save(y_qualtrics_dat_sel, file = paste0(selected_clean_dat_path, "y_qualtrics_dat_sel.RData"))
save(p_qualtrics_dat_sel, file = paste0(selected_clean_dat_path, "p_qualtrics_dat_sel.RData"))