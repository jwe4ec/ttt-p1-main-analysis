# ---------------------------------------------------------------------------- #
# Extract Person-Specific Temporal Network Results -----
# Authors: Josip Razum, Jeremy W. Eberle
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

# Load package

groundhog.library("dplyr", groundhog_day)

# ---------------------------------------------------------------------------- #
# Import results ----
# ---------------------------------------------------------------------------- #

varfit           <- readRDS(file = "./02_networks/results/raw/varfit.RDS")
varfit_control   <- readRDS(file = "./02_networks/results/raw/varfit_control.RDS")
varfit_fun       <- readRDS(file = "./02_networks/results/raw/varfit_fun.RDS")

mlvarfit         <- readRDS(file = "./02_networks/results/raw/mlvarfit.RDS")
mlvarfit_control <- readRDS(file = "./02_networks/results/raw/mlvarfit_control.RDS")
mlvarfit_fun     <- readRDS(file = "./02_networks/results/raw/mlvarfit_fun.RDS")

# Load Lifepak IDs

load("./02_networks/data/final_clean/ids.RDS")

# ---------------------------------------------------------------------------- #
# Label idiographic VAR fit list elements by "lifepak_id" ----
# ---------------------------------------------------------------------------- #

# Define function to label idiographic VAR fit list elements by "lifepak_id"

name_var_fit_elements <- function(fit, lifepak_ids) {
  if (length(fit) == length(lifepak_ids)) {
    names(fit) <- lifepak_ids
  } else {
    stop("Length of VAR fit list differs from length of Lifepak IDs")
  }
  
  return(fit)
}

# Run function for VAR fit list elements

varfit         <- name_var_fit_elements(varfit,         ids)
varfit_control <- name_var_fit_elements(varfit_control, ids)
varfit_fun     <- name_var_fit_elements(varfit_fun,     ids)

# ---------------------------------------------------------------------------- #
# Check that ML-VAR clusters are Lifepak IDs ----
# ---------------------------------------------------------------------------- #

# All clusters are Lifepak IDs but in a different order (will reorder later in script)

length(intersect(ids, unique(mlvarfit$parameters$wilevel.standardized$stdyx.standardized$cluster))) == length(ids)

ids
unique(mlvarfit$parameters$wilevel.standardized$stdyx.standardized$cluster)

# ---------------------------------------------------------------------------- #
# Define parameter numbers for temporal results for VAR models ----
# ---------------------------------------------------------------------------- #

# Manually define parameter numbers in fit object for standardized autoregressive 
# and cross-lagged model coefficients

  # For 8-node model (parameters 1-64, whereas 65-92 are contemporaneous effects, 
  # 93-100 are intercepts, and 101-108 are residual variances)

# View(varfit[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var <- 1:64

  # For 7-node model with "control" (parameters 1-49, whereas 50-70 are contemporaneous 
  # effects, 71-77 are intercepts, and 78-84 are residual variances)

# View(varfit_control[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var_control <- 1:49

  # For 7-node model with "fun" (same as for 7-node model with "control")

# View(varfit_fun[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var_fun <- 1:49

# ---------------------------------------------------------------------------- #
# Extract and threshold temporal results for VAR models ----
# ---------------------------------------------------------------------------- #

# Define function to extract and threshold standardized temporal results

extract_temporal_var_results <- function(fit, param_numbers) {
  # Extract standardized autoregressive and cross-lagged model coefficients
  
  results <- lapply(fit, function(x) {
    out <- cbind(param_numbers,
                 x$parameters$stdyx.standardized[param_numbers,
                                                 c("paramHeader", "param", 
                                                   "est", "lower_2.5ci", "upper_2.5ci")])
    
    names(out)[names(out) == "param_numbers"] <- "param_number"
    
    return(out)
  })
  
  # Add "lifepak_id" column to each data frame in list
  
  for (i in 1:length(results)) {
    results[[i]]$lifepak_id <- names(results[i])
    
    results[[i]] <- results[[i]][, c("lifepak_id", names(results[[i]])[names(results[[i]]) != "lifepak_id"])]
  }
  
  # Stack results for all participants in one data frame
  
  results <- do.call(rbind, results)
  
  # Create columns for criterion and predictor variable names (based on all 8
  # possible nodes; in 7-node networks, only 7 of the 8 lines below will apply)
  
  results$criterion <- NA
  results$criterion[results$paramHeader == "BAD.ON"]  <- "bad"
  results$criterion[results$paramHeader == "CONT.ON"] <- "control"
  results$criterion[results$paramHeader == "ENER.ON"] <- "energy"
  results$criterion[results$paramHeader == "FOC.ON"]  <- "focus"
  results$criterion[results$paramHeader == "FUN.ON"]  <- "fun"
  results$criterion[results$paramHeader == "INT.ON"]  <- "interest"
  results$criterion[results$paramHeader == "MOVE.ON"] <- "movement"
  results$criterion[results$paramHeader == "SAD.ON"]  <- "sad"
  
  results$predictor <- NA
  results$predictor[results$param == "BAD&1"]  <- "bad"
  results$predictor[results$param == "CONT&1"] <- "control"
  results$predictor[results$param == "ENER&1"] <- "energy"
  results$predictor[results$param == "FOC&1"]  <- "focus"
  results$predictor[results$param == "FUN&1"]  <- "fun"
  results$predictor[results$param == "INT&1"]  <- "interest"
  results$predictor[results$param == "MOVE&1"] <- "movement"
  results$predictor[results$param == "SAD&1"]  <- "sad"
  
  # Threshold edges based on significance
  
  results$est_thres <- NA
  
  results$est_thres <- ifelse((results$lower_2.5ci > 0 & results$upper_2.5ci > 0) |
                                (results$lower_2.5ci < 0 & results$upper_2.5ci < 0), results$est, 0)
  
  return(results)
}

# Run function

results_var         <- extract_temporal_var_results(varfit,         param_numbers_var)
results_var_control <- extract_temporal_var_results(varfit_control, param_numbers_var_control)
results_var_fun     <- extract_temporal_var_results(varfit_fun,     param_numbers_var_fun)

# ---------------------------------------------------------------------------- #
# Define parameter numbers for temporal results for ML-VAR models ----
# ---------------------------------------------------------------------------- #

# Manually define parameter numbers in fit object for standardized autoregressive 
# and cross-lagged model coefficients

  # For 8-node model (parameters 1-64, whereas 65-92 are contemporaneous effects 
  # and 93-100 are residual variances)

# View(mlvarfit$parameters$wilevel.standardized$stdyx.standardized)
param_numbers_mlvar <- 1:64

  # For 7-node model with "control" (parameters 1-49, whereas 50-70 are contemporaneous 
  # effects and 71-77 are residual variances)

# View(mlvarfit_control$parameters$wilevel.standardized$stdyx.standardized)
param_numbers_mlvar_control <- 1:49

  # For 7-node model with "fun" (same as for 7-node model with "control")

# View(mlvarfit_fun$parameters$wilevel.standardized$stdyx.standardized)
param_numbers_mlvar_fun <- 1:49

# ---------------------------------------------------------------------------- #
# Extract and threshold temporal results for ML-VAR models ----
# ---------------------------------------------------------------------------- #

# Define function to extract and threshold standardized temporal results

extract_temporal_mlvar_results <- function(fit, param_numbers) {
  # Extract standardized autoregressive and cross-lagged model coefficients
  
  results <- fit$parameters$wilevel.standardized$stdyx.standardized
  
  names(results)[names(results) == "cluster"] <- "lifepak_id"
  
  # Note: Use factor version of "lifepak_id" to preserve participant order
  
  results$lifepak_id_fct <- factor(results$lifepak_id, levels = unique(results$lifepak_id))
  
  results <- results %>%
    group_by(lifepak_id_fct) %>%
    slice(param_numbers) %>%
    ungroup()
  
  results <- results[, c("lifepak_id", "paramHeader", "param", 
                         "est", "lower_2.5ci", "upper_2.5ci")]
  
  # Create columns for criterion and predictor variable names (based on all 8
  # possible nodes; in 7-node networks, only 7 of the 8 lines below will apply)
  
  results$criterion <- NA
  results$criterion[grepl("BAD.ON",  results$paramHeader)] <- "bad"
  results$criterion[grepl("CONT.ON", results$paramHeader)] <- "control"
  results$criterion[grepl("ENER.ON", results$paramHeader)] <- "energy"
  results$criterion[grepl("FOC.ON",  results$paramHeader)] <- "focus"
  results$criterion[grepl("FUN.ON",  results$paramHeader)] <- "fun"
  results$criterion[grepl("INT.ON",  results$paramHeader)] <- "interest"
  results$criterion[grepl("MOVE.ON", results$paramHeader)] <- "movement"
  results$criterion[grepl("SAD.ON",  results$paramHeader)] <- "sad"
  
  results$predictor <- NA
  results$predictor[results$param == "BAD&1"]  <- "bad"
  results$predictor[results$param == "CONT&1"] <- "control"
  results$predictor[results$param == "ENER&1"] <- "energy"
  results$predictor[results$param == "FOC&1"]  <- "focus"
  results$predictor[results$param == "FUN&1"]  <- "fun"
  results$predictor[results$param == "INT&1"]  <- "interest"
  results$predictor[results$param == "MOVE&1"] <- "movement"
  results$predictor[results$param == "SAD&1"]  <- "sad"
  
  # Threshold edges based on significance
  
  results$est_thres <- NA
  
  results$est_thres <- ifelse((results$lower_2.5ci > 0 & results$upper_2.5ci > 0) |
                                (results$lower_2.5ci < 0 & results$upper_2.5ci < 0), results$est, 0)

  return(results)
}

# Run function

results_mlvar         <- extract_temporal_mlvar_results(mlvarfit,         param_numbers_mlvar)
results_mlvar_control <- extract_temporal_mlvar_results(mlvarfit_control, param_numbers_mlvar_control)
results_mlvar_fun     <- extract_temporal_mlvar_results(mlvarfit_fun,     param_numbers_mlvar_fun)

# ---------------------------------------------------------------------------- #
# Create thresholded adjacency matrices for VAR and ML-VAR models ----
# ---------------------------------------------------------------------------- #

# Define function to create directed adjacency matrix of thresholded autoregressive 
# and cross-lagged coefficients for a given participant (where rows are predictors
# and columns are criterions)

create_thres_adj_mat <- function(participant_results) {
  variables <- unique(c(participant_results$criterion, participant_results$predictor))
  
  adj_mat <- matrix(NA, nrow = length(variables), ncol = length(variables))
  
  rownames(adj_mat) <- variables
  colnames(adj_mat) <- variables
  
  for (i in 1:nrow(participant_results)) {
    predictor <- participant_results$predictor[i]
    criterion <- participant_results$criterion[i]
    relationship <- participant_results$est_thres[i]
    adj_mat[predictor, criterion] <- relationship
  }
  
  return(adj_mat)
}

# Define function to group results by participant and create adjacency matrices

create_thres_adj_mats <- function(results) {
  # Note: Use factor version of "lifepak_id" to preserve participant order
  
  results$lifepak_id_fct <- factor(results$lifepak_id, levels = unique(results$lifepak_id))
  
  thres_adj_mats <- results %>%
    split(.$lifepak_id_fct) %>%
    lapply(create_thres_adj_mat)
  
  return(thres_adj_mats)
}

# Run function to create adjacency matrices

thres_adj_mats_var           <- create_thres_adj_mats(results_var)
thres_adj_mats_var_control   <- create_thres_adj_mats(results_var_control)
thres_adj_mats_var_fun       <- create_thres_adj_mats(results_var_fun)

thres_adj_mats_mlvar         <- create_thres_adj_mats(results_mlvar)
thres_adj_mats_mlvar_control <- create_thres_adj_mats(results_mlvar_control)
thres_adj_mats_mlvar_fun     <- create_thres_adj_mats(results_mlvar_fun)

# ---------------------------------------------------------------------------- #
# Sort results and adjacency matrices by "lifepak_id" ----
# ---------------------------------------------------------------------------- #

# VAR results and adjacency matrices are already sorted by "lifepak_id"

all(unique(results_var$lifepak_id)         == ids)
all(unique(results_var_control$lifepak_id) == ids)
all(unique(results_var_fun$lifepak_id)     == ids)

all(names(thres_adj_mats_var)         == ids)
all(names(thres_adj_mats_var_control) == ids)
all(names(thres_adj_mats_var_fun)     == ids)

# Sort ML-VAR results and adjacency matrices

results_mlvar         <- results_mlvar[order(results_mlvar$lifepak_id), ]
results_mlvar_control <- results_mlvar_control[order(results_mlvar_control$lifepak_id), ]
results_mlvar_fun     <- results_mlvar_fun[order(results_mlvar_fun$lifepak_id), ]

thres_adj_mats_mlvar         <- thres_adj_mats_mlvar[as.character(sort(as.integer(names(thres_adj_mats_mlvar))))]
thres_adj_mats_mlvar_control <- thres_adj_mats_mlvar_control[as.character(sort(as.integer(names(thres_adj_mats_mlvar_control))))]
thres_adj_mats_mlvar_fun     <- thres_adj_mats_mlvar_fun[as.character(sort(as.integer(names(thres_adj_mats_mlvar_fun))))]

# ---------------------------------------------------------------------------- #
# Export temporal network results and adjacency matrices ----
# ---------------------------------------------------------------------------- #

extracted_results_path <- "./02_networks/results/extracted/"

dir.create(extracted_results_path)

save(results_var,           file = paste0(extracted_results_path, "results_var.RDS"))
save(results_var_control,   file = paste0(extracted_results_path, "results_var_control.RDS"))
save(results_var_fun,       file = paste0(extracted_results_path, "results_var_fun.RDS"))

save(results_mlvar,         file = paste0(extracted_results_path, "results_mlvar.RDS"))
save(results_mlvar_control, file = paste0(extracted_results_path, "results_mlvar_control.RDS"))
save(results_mlvar_fun,     file = paste0(extracted_results_path, "results_mlvar_fun.RDS"))

save(thres_adj_mats_var,           file = paste0(extracted_results_path, "thres_adj_mats_var.RDS"))
save(thres_adj_mats_var_control,   file = paste0(extracted_results_path, "thres_adj_mats_var_control.RDS"))
save(thres_adj_mats_var_fun,       file = paste0(extracted_results_path, "thres_adj_mats_var_fun.RDS"))

save(thres_adj_mats_mlvar,         file = paste0(extracted_results_path, "thres_adj_mats_mlvar.RDS"))
save(thres_adj_mats_mlvar_control, file = paste0(extracted_results_path, "thres_adj_mats_mlvar_control.RDS"))
save(thres_adj_mats_mlvar_fun,     file = paste0(extracted_results_path, "thres_adj_mats_mlvar_fun.RDS"))