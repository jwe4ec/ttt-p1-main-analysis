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
# Compute raw change scores for Qualtrics outcomes ----
# ---------------------------------------------------------------------------- #

comb_dat$y_cdi_m_chg     <- comb_dat$y3m_cdi_mean     - comb_dat$yb_cdi_mean
comb_dat$p_cdi_m_chg     <- comb_dat$p3m_cdi_mean     - comb_dat$pb_cdi_mean
comb_dat$y_bhs_m_chg     <- comb_dat$y3m_bhs_mean     - comb_dat$yb_bhs_mean
comb_dat$y_pcsc_m_chg    <- comb_dat$y3m_pcsc_mean    - comb_dat$yb_pcsc_mean
comb_dat$y_bads_ac_m_chg <- comb_dat$y3m_bads_ac_mean - comb_dat$yb_bads_ac_mean
comb_dat$y_bads_ar_m_chg <- comb_dat$y3m_bads_ar_mean - comb_dat$yb_bads_ar_mean

# ---------------------------------------------------------------------------- #
# TODO: Restrict to participants with complete Qualtrics outcome data ----
# ---------------------------------------------------------------------------- #

# TODO: Temporarily restrict to participants with complete Qualtrics outcome data
# at baseline and 3 months (once networks are refit, this will already have been done
# and "comb_dat" will have the final analysis sample based on "net_params_var_mlvar")

load(paste0("./02_networks/data/merged_clean/lifepak_ids_qualtrics_compl.RData"))

comb_dat <- comb_dat[comb_dat$lifepak_id %in% lifepak_ids_qualtrics_compl, ]





# ---------------------------------------------------------------------------- #
# Identify names of predictors and outcomes of prediction models ----
# ---------------------------------------------------------------------------- #

# Identify names of raw change scores

m_chg_outcomes <- c("y_cdi_m_chg", "p_cdi_m_chg", "y_bhs_m_chg", "y_pcsc_m_chg", 
                    "y_bads_ac_m_chg", "y_bads_ar_m_chg")

length(m_chg_outcomes) == 6

all(m_chg_outcomes %in% names(comb_dat))

# Identify names of Qualtrics outcomes at baseline

b_outcomes <- c("yb_cdi_mean", "pb_cdi_mean", "yb_bhs_mean", "yb_pcsc_mean", 
                "yb_bads_ac_mean", "yb_bads_ar_mean")

length(b_outcomes) == 6

all(b_outcomes %in% names(qualtrics_dat_items_scales_dem))

# Identify names of raw means of EMA items

raw_means_ema_items <- paste0(node_vars, "_m")

length(raw_means_ema_items) == 8

setequal(c("lifepak_id", raw_means_ema_items), names(raw_means))

# Identify names of centrality and density parameters from VAR and ML-VAR models

density_params <- c("inter_conn__var",   "intra_conn__var",
                    "inter_conn__mlvar", "intra_conn__mlvar")

density_params_var   <- density_params[grepl("__var",   density_params)]
density_params_mlvar <- density_params[grepl("__mlvar", density_params)]

length(density_params)       == 4
length(density_params_var)   == 2
length(density_params_mlvar) == 2

centrality_params <- setdiff(names(net_params_var_mlvar), c("lifepak_id", density_params))

centrality_params_var   <- centrality_params[grepl("__var",   centrality_params)]
centrality_params_mlvar <- centrality_params[grepl("__mlvar", centrality_params)]

length(centrality_params)       == 88
length(centrality_params_var)   == 44
length(centrality_params_mlvar) == 44

all_net_params <- c(density_params_var, density_params_mlvar, centrality_params_var, centrality_params_mlvar)

setequal(c("lifepak_id", all_net_params), names(net_params_var_mlvar))

# ---------------------------------------------------------------------------- #
# Select relevant columns for prediction models ----
# ---------------------------------------------------------------------------- #

# Select only "lifepak_id", predictors, and outcomes

pred_dat <- comb_dat[, c("lifepak_id", m_chg_outcomes, b_outcomes, raw_means_ema_items, all_net_params)]

# ---------------------------------------------------------------------------- #
# TODO: Select relevant columns for demographics table ----
# ---------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
# TODO: Select relevant columns for missing data rates ----
# ---------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
# Save data ----
# ---------------------------------------------------------------------------- #

final_clean_path <- "./02_networks/data/final_clean/"

save(comb_dat, file = paste0(final_clean_path, "comb_dat.RData"))
save(pred_dat, file = paste0(final_clean_path, "pred_dat.RData"))