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

# Load packages

groundhog.library("stringr", groundhog_day)

# ---------------------------------------------------------------------------- #
# Import clean EMA data, computed network parameters, and merged Qualtrics data  ----
# ---------------------------------------------------------------------------- #

# TODO: Update data after Phase I data cleaning is complete





load("./02_networks/data/final_clean/data_var_qualtrics_compl.RDS")
ema_dat <- data_var_qualtrics_compl

# TODO: Refit network models after Phase I data cleaning is complete





load("./02_networks/results/net_params/net_params_var_mlvar.RDS")

# TODO: Import merged Qualtrics data (update once data cleaning is complete)

load("./02_networks/data/merged_clean/qualtrics_dat_items_scales_dem.RData")





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

# TODO: Temporarily recode 5-digit LifePak IDs in "net_params_var_mlvar" to contain 
# leading 0 so that all LifePak IDs are 6 digits (this is Isaac's format)

net_params_var_mlvar$lifepak_id <- str_pad(net_params_var_mlvar$lifepak_id, width = 6, side = "left", pad = "0")





# Add raw means of EMA items to computed network parameters

comb_dat <- merge(net_params_var_mlvar, raw_means,          "lifepak_id", all.x = TRUE, sort = FALSE)

# Add merged Qualtrics data

comb_dat <- merge(comb_dat, qualtrics_dat_items_scales_dem, "lifepak_id", all.x = TRUE, sort = FALSE)

# ---------------------------------------------------------------------------- #
# TODO: Restrict to participants with complete Qualtrics outcome data ----
# ---------------------------------------------------------------------------- #

# TODO: Temporarily restrict to participants with complete Qualtrics outcome data
# at baseline and 3 months (once networks are refit, this will already have been done
# and "comb_dat" will have the final analysis sample based on "net_params_var_mlvar")

load(paste0("./02_networks/data/merged_clean/lifepak_ids_qualtrics_compl.RData"))

comb_dat <- comb_dat[comb_dat$lifepak_id %in% lifepak_ids_qualtrics_compl, ]





# ---------------------------------------------------------------------------- #
# TODO: Select relevant columns for prediction models ----
# ---------------------------------------------------------------------------- #

# TODO: Restrict to only "lifepak_id", predictors, and outcomes





# ---------------------------------------------------------------------------- #
# TODO: Select relevant columns for demographics table ----
# ---------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
# TODO: Select relevant columns for missing data rates ----
# ---------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
# TODO: Save data ----
# ---------------------------------------------------------------------------- #




