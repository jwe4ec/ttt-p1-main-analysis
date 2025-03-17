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
# Merge datasets ----
# ---------------------------------------------------------------------------- #

# TODO




