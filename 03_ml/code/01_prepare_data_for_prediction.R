# ---------------------------------------------------------------------------- #
# Prepare Data for Prediction Models -----
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
# Import clean EMA data, computed network parameters, and clean Qualtrics data  ----
# ---------------------------------------------------------------------------- #

# TODO: Update data after Phase I data cleaning is complete





load("./02_networks/data/final_clean/data_var.RDS")
ema_dat <- data_var

# TODO: Refit network models after Phase I data cleaning is complete





load("./02_networks/results/net_params/net_params_var_mlvar.RDS")

# TODO: Finalize file organization after Phase I data cleaning is complete





clean_qualtrics_dat_path <- "R:/MSS/Schleider_Lab/jslab/TRACK to TREAT/Data/Clean Data (Isaac)/"

y_qualtrics_dat <- readRDS(paste0(clean_qualtrics_dat_path, "Phase 1 Youth Qualtrics Data.RDS"))
p_qualtrics_dat <- readRDS(paste0(clean_qualtrics_dat_path, "Phase 1 Parent Qualtrics Data.RDS"))

# ---------------------------------------------------------------------------- #
# Compute raw means of clean EMA items across time per participant ----
# ---------------------------------------------------------------------------- #

node_vars <- c("bad", "control", "energy", "focus", "fun", "interest", "movement", "sad")

raw_means <- data.frame(lifepak_id = unique(ema_dat$lifepak_id))

for (node_var in node_vars) {
  node_var_m <- paste0(node_var, "_m")
  
  ag <- aggregate(ema_dat[, node_var], list(lifepak_id = ema_dat$lifepak_id), mean, na.rm = TRUE)
  names(ag)[names(ag) == "x"] <- node_var_m

  raw_means <- merge(raw_means, ag, "lifepak_id", all.x = TRUE, sort = FALSE)
}

# ---------------------------------------------------------------------------- #
# Define and restrict to Qualtrics columns of interest ----
# ---------------------------------------------------------------------------- #

# Define metadata columns

y_meta_cols <- c("lsmh_id", "lifepak_id", "yb_administration", "yb_date", "y3m_date")

p_meta_cols <- c("lsmh_id", "pb_administration", "pb_date", "p3m_date")

all(y_meta_cols %in% names(y_qualtrics_dat))
all(p_meta_cols %in% names(p_qualtrics_dat))

# TODO (Check factor structure): Define relevant outcome columns at baseline and 3 months. 
# Note: Reverse-coded columns have already been unreversed (as part of data cleaning).





  # Youth-rated

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

all(c(y_cdi_cols, y_bhs_cols, y_pcsc_cols, y_bads_ac_cols, y_bads_ar_cols) %in% names(y_qualtrics_dat))

# Parent-rated

p_cdi_cols     <- c(paste0("pb_cdi_",  1:17),
                    paste0("p3m_cdi_", 1:17),
                    "pb_cdi_mean", "p3m_cdi_mean")

  # TODO: Add parent demographics once cleaned (at minimum, parent age, sex, gender, 
  # race [if assessed], ethnicity, single- vs. co-parent, depression symptom severity)





p_dem_cols     <- c("pb_childage", "pb_childsex", "pb_childgender", "pb_childethnicity", 
                    "pb_birthorder", "pb_n_sisters", "pb_n_brothers", "pb_grade", 
                    "pb_school", "pb_income", "pb_dependent", "pb_childtx_lifetime", 
                    "p3m_childtx_lifetime", "pb_childtx_current")

all(c(p_cdi_cols, p_dem_cols) %in% names(p_qualtrics_dat))

# TODO: Restrict to columns of interest (separate outcomes from demographics)





# ---------------------------------------------------------------------------- #
# TODO: Merge datasets ----
# ---------------------------------------------------------------------------- #

comb_dat <- merge(net_params_var_mlvar, raw_means, "lifepak_id", all.x = TRUE, sort = FALSE)




# ---------------------------------------------------------------------------- #
# TODO: Restrict to participants with complete baseline and 3-month data ----
# ---------------------------------------------------------------------------- #




