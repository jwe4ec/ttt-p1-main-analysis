# ---------------------------------------------------------------------------- #
# Further Clean Qualtrics Data -----
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
# Import selected clean Qualtrics data  ----
# ---------------------------------------------------------------------------- #

# TODO: Finalize file organization after Phase I data cleaning is complete






selected_clean_dat_path <- "./02_networks/data/selected_clean/"

load(paste0(selected_clean_dat_path, "y_qualtrics_dat_sel.RData"))
load(paste0(selected_clean_dat_path, "p_qualtrics_dat_sel.RData"))

# ---------------------------------------------------------------------------- #
# Merge youth and parent Qualtrics datasets ----
# ---------------------------------------------------------------------------- #

all(y_qualtrics_dat_sel$lsmh_id == p_qualtrics_dat_sel$lsmh_id)

qualtrics_dat_items_scales_dem <- merge(y_qualtrics_dat_sel, p_qualtrics_dat_sel, "lsmh_id", all.x = TRUE)

# 106 participants enrolled

length(qualtrics_dat_items_scales_dem$lifepak_id) == 106

# ---------------------------------------------------------------------------- #
# Identify participants with complete Qualtrics outcome data ----
# ---------------------------------------------------------------------------- #

# Identify participants with complete data on average item scores for all relevant 
# outcomes at baseline and 3 months

mean_cols_b  <- c("yb_cdi_mean",  "yb_bhs_mean",  "yb_pcsc_mean",  "yb_bads_ac_mean", 
                  "yb_bads_ar_mean",  "pb_cdi_mean")
mean_cols_3m <- c("y3m_cdi_mean", "y3m_bhs_mean", "y3m_pcsc_mean", "y3m_bads_ac_mean",
                  "y3m_bads_ar_mean", "p3m_cdi_mean")

mean_cols <- c(mean_cols_b, mean_cols_3m)

# TODO (temporarily do this; asked Isaac to fold it into data cleaning): Change
# NaN for mean scores to NA





sum(is.na(qualtrics_dat_items_scales_dem[, mean_cols]))            == 94
sum(apply(qualtrics_dat_items_scales_dem[, mean_cols], 2, is.nan)) == 94

qualtrics_dat_items_scales_dem[, mean_cols][is.na(qualtrics_dat_items_scales_dem[, mean_cols])] <- NA

sum(is.na(qualtrics_dat_items_scales_dem[, mean_cols]))            == 94
sum(apply(qualtrics_dat_items_scales_dem[, mean_cols], 2, is.nan)) == 0

# Compute indicator of complete data for mean columns at baseline and 3 months

qualtrics_dat_items_scales_dem$mean_cols_b_compl <- NA
qualtrics_dat_items_scales_dem$mean_cols_3m_compl <- NA

qualtrics_dat_items_scales_dem$mean_cols_b_compl <-
  apply(qualtrics_dat_items_scales_dem[, mean_cols_b],  1, function(x) as.integer(all(!is.na(x))))
qualtrics_dat_items_scales_dem$mean_cols_3m_compl <-
  apply(qualtrics_dat_items_scales_dem[, mean_cols_3m], 1, function(x) as.integer(all(!is.na(x))))

# 1 participant is missing data for 1 mean score at baseline

sum(qualtrics_dat_items_scales_dem$mean_cols_b_compl == 0) == 1
sum(is.na(qualtrics_dat_items_scales_dem[, mean_cols_b]))  == 1

# 21 participants are missing data for 93 mean scores at 3 months

sum(qualtrics_dat_items_scales_dem$mean_cols_3m_compl == 0) == 21
sum(is.na(qualtrics_dat_items_scales_dem[, mean_cols_3m]))  == 93

# 22 participants are missing data for at least 1 mean score at baseline or 3 months,
# leaving 84 participants with complete outcome data

lifepak_ids_qualtrics_incompl <- 
  qualtrics_dat_items_scales_dem$lifepak_id[qualtrics_dat_items_scales_dem$mean_cols_b_compl == 0 |
                                              qualtrics_dat_items_scales_dem$mean_cols_3m_compl == 0]

length(lifepak_ids_qualtrics_incompl) == 22

lifepak_ids_qualtrics_compl <- setdiff(qualtrics_dat_items_scales_dem$lifepak_id, 
                                       lifepak_ids_qualtrics_incompl)

length(lifepak_ids_qualtrics_compl) == 84

# ---------------------------------------------------------------------------- #
# Save data and IDs  ----
# ---------------------------------------------------------------------------- #

# TODO: Save merged Qualtrics data and LifePak IDs of participants who have complete
# Qualtrics outcome data




