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

# Load packages

pkgs <- c("dplyr", "tidyr", "Hmisc", "networktools", "netcontrol") # TODO: Determine what packages are needed

groundhog.library(pkgs, groundhog_day)

# ---------------------------------------------------------------------------- #
# Import results ----
# ---------------------------------------------------------------------------- #

varfit           <- readRDS(file = "./02_networks/results/varfit.RDS")
varfit_control   <- readRDS(file = "./02_networks/results/varfit_control.RDS")
varfit_fun       <- readRDS(file = "./02_networks/results/varfit_fun.RDS")

mlvarfit         <- readRDS(file = "./02_networks/results/mlvarfit.RDS")
mlvarfit_control <- readRDS(file = "./02_networks/results/mlvarfit_control.RDS")
mlvarfit_fun     <- readRDS(file = "./02_networks/results/mlvarfit_fun.RDS")

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
# TODO: Check that ML-VAR clusters correspond to Lifepak IDs ----
# ---------------------------------------------------------------------------- #

# TODO: Although all ML-VAR clusters reflect Lifepak IDs, the order differs

length(intersect(ids, unique(mlvarfit$parameters$wilevel.standardized$stdyx.standardized$cluster))) == length(ids)

ids
unique(mlvarfit$parameters$wilevel.standardized$stdyx.standardized$cluster)





# ---------------------------------------------------------------------------- #
# Define parameter numbers for temporal results for VAR models ----
# ---------------------------------------------------------------------------- #

# Manually define parameter numbers in fit object for standardized autoregressive 
# and cross-lagged model coefficients

  # For 8-node VAR model (parameters 1-64, whereas 65-92 are contemporaneous effects, 
  # 93-100 are intercepts, and 101-108 are residual variances)

# View(varfit[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var <- 1:64

  # For 7-node VAR model with "control" (parameters 1-49, whereas 50-70 are contemporaneous 
  # effects, 71-77 are intercepts, and 78-84 are residual variances)

# View(varfit_control[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var_control <- 1:49

  # For 7-node VAR model with "fun" (same as for 7-node VAR model with "control")

# View(varfit_fun[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var_fun <- 1:49

# ---------------------------------------------------------------------------- #
# Extract and threshold temporal results for VAR models ----
# ---------------------------------------------------------------------------- #

# Define function to extract and threshold standardized temporal results

extract_temporal_results <- function(fit, param_numbers) {
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

results_var         <- extract_temporal_results(varfit,         param_numbers_var)
results_var_control <- extract_temporal_results(varfit_control, param_numbers_var_control)
results_var_fun     <- extract_temporal_results(varfit_fun,     param_numbers_var_fun)

# ---------------------------------------------------------------------------- #
# Define parameter numbers for temporal results for VAR models ----
# ---------------------------------------------------------------------------- #

# Manually define parameter numbers in fit object for standardized autoregressive 
# and cross-lagged model coefficients

  # For 8-node VAR model (parameters 1-64, whereas 65-92 are contemporaneous effects, 
  # 93-100 are intercepts, and 101-108 are residual variances)

# View(varfit[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var <- 1:64

  # For 7-node VAR model with "control" (parameters 1-49, whereas 50-70 are contemporaneous 
  # effects, 71-77 are intercepts, and 78-84 are residual variances)

# View(varfit_control[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var_control <- 1:49

  # For 7-node VAR model with "fun" (same as for 7-node VAR model with "control")

# View(varfit_fun[[1]]$parameters$stdyx.standardized) # Example participant
param_numbers_var_fun <- 1:49

# ---------------------------------------------------------------------------- #
# Extract and threshold temporal results for VAR models ----
# ---------------------------------------------------------------------------- #

# Define function to extract and threshold standardized temporal results

extract_temporal_results <- function(fit, param_numbers) {
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

results_var         <- extract_temporal_results(varfit,         param_numbers_var)
results_var_control <- extract_temporal_results(varfit_control, param_numbers_var_control)
results_var_fun     <- extract_temporal_results(varfit_fun,     param_numbers_var_fun)

# ---------------------------------------------------------------------------- #
# Create thresholded adjacency matrices for VAR models ----
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

thres_adj_mats_var         <- create_thres_adj_mats(results_var)
thres_adj_mats_var_control <- create_thres_adj_mats(results_var_control)
thres_adj_mats_var_fun     <- create_thres_adj_mats(results_var_fun)




# TODO: Separate extraction of model coefficients and creation of adjacency
# matrices from computation of network parameters. Also, condense code using
# custom functions rather than repeating code across several models.




# ---------------------------------------------------------------------------- #
# Extract and threshold temporal results for 8-node ML-VAR model ----
# ---------------------------------------------------------------------------- #

param_numbers_mlvar <- 1:64

# Extract standardized autoregressive and cross-lagged model coefficients

results_mlvar <- mlvarfit$parameters$wilevel.standardized$stdyx.standardized

names(results_mlvar)[names(results_mlvar) == "cluster"] <- "lifepak_id"

results_mlvar <- results_mlvar %>%
  group_by(lifepak_id) %>%
  slice(param_numbers_mlvar) %>%
  ungroup()

results_mlvar <- results_mlvar[, c("lifepak_id", "paramHeader", "param", 
                                   "est", "lower_2.5ci", "upper_2.5ci")]

# Create columns for criterion and predictor variable names (based on all 8
# possible nodes; in 7-node networks, only 7 of the 8 lines below will apply)

results_mlvar$criterion <- NA
results_mlvar$criterion[grepl("BAD.ON",  results_mlvar$paramHeader)] <- "bad"
results_mlvar$criterion[grepl("CONT.ON", results_mlvar$paramHeader)] <- "control"
results_mlvar$criterion[grepl("ENER.ON", results_mlvar$paramHeader)] <- "energy"
results_mlvar$criterion[grepl("FOC.ON",  results_mlvar$paramHeader)] <- "focus"
results_mlvar$criterion[grepl("FUN.ON",  results_mlvar$paramHeader)] <- "fun"
results_mlvar$criterion[grepl("INT.ON",  results_mlvar$paramHeader)] <- "interest"
results_mlvar$criterion[grepl("MOVE.ON", results_mlvar$paramHeader)] <- "movement"
results_mlvar$criterion[grepl("SAD.ON",  results_mlvar$paramHeader)] <- "sad"

results_mlvar$predictor <- NA
results_mlvar$predictor[results_mlvar$param == "BAD&1"]  <- "bad"
results_mlvar$predictor[results_mlvar$param == "CONT&1"] <- "control"
results_mlvar$predictor[results_mlvar$param == "ENER&1"] <- "energy"
results_mlvar$predictor[results_mlvar$param == "FOC&1"]  <- "focus"
results_mlvar$predictor[results_mlvar$param == "FUN&1"]  <- "fun"
results_mlvar$predictor[results_mlvar$param == "INT&1"]  <- "interest"
results_mlvar$predictor[results_mlvar$param == "MOVE&1"] <- "movement"
results_mlvar$predictor[results_mlvar$param == "SAD&1"]  <- "sad"

# Threshold edges based on significance

results_mlvar$est_thres <- NA

results_mlvar$est_thres <- ifelse((results_mlvar$lower_2.5ci > 0 & results_mlvar$upper_2.5ci > 0) |
                                    (results_mlvar$lower_2.5ci < 0 & results_mlvar$upper_2.5ci < 0), results_mlvar$est, 0)

# TODO: Create adjacency matrices using functions for VAR models above





# ---------------------------------------------------------------------------- #
# Extracting the person-specific parameters for the 7-node ML-VAR network with control ----
# ---------------------------------------------------------------------------- #

parameters_mlvar_solution_control <- mlvarfit_control$parameters$wilevel.standardized$stdyx.standardized

results_mlvar_control <- parameters_mlvar_solution_control

# Subset the dataframe to keep only the desired rows
results_mlvar_control <- results_mlvar_control %>%
  group_by(cluster) %>%
  slice_head(n = 49) %>%
  ungroup()

#create a column in the dataframe that signifies to which participant do the given set of coefficients belong. 
# Calculate the number of participants
num_participants <- ceiling(nrow(results_mlvar_control) / 49)

# Create the participant labels
participant_labels <- rep(paste0("Participant ", 1:num_participants), each = 49, length.out = nrow(results_mlvar_control))

# Add the new column to the dataframe
results_mlvar_control$Participant <- participant_labels

results_mlvar_control <- results_mlvar_control[, c("Participant", "paramHeader", "param", "est", "lower_2.5ci", "upper_2.5ci", "sig", "cluster")]

#create the column with criterion variable names

variables <- c("bad", "control", "energy", "focus", "interest", "movement", "sad")

repeated_variables <- rep(rep(variables, each = 7), length.out = nrow(results_mlvar_control))

results_mlvar_control$criterion <- repeated_variables

#create the column with predictor variable names
repeated_predictor_variables <- rep(variables, length.out = nrow(results_mlvar_control))

results_mlvar_control$predictor <- repeated_predictor_variables

#replacing non-significant edges with zeroes

results_mlvar_control$est <- ifelse((results_mlvar_control$lower_2.5ci > 0 & results_mlvar_control$upper_2.5ci > 0) |
                              (results_mlvar_control$lower_2.5ci < 0 & results_mlvar_control$upper_2.5ci < 0), results_mlvar_control$est, 0)

# ---------------------------------------------------------------------------- #
# Creating the adjacency matrix from ML-VAR parameters for the 7-node ML-VAR network with control ----
# ---------------------------------------------------------------------------- #

#create an adjacency matrix of auto-regressive and cross-regressive coefficients
#create a directed adjacency matrix for each participant

# Function to create adjacency matrix for a given participant's data
create_adjacency_matrix <- function(participant_data) {
  variables_mlvar_control <- unique(c(participant_data$criterion, participant_data$predictor))
  adj_matrix <- matrix(NA, nrow = length(variables_mlvar_control), ncol = length(variables_mlvar_control))
  rownames(adj_matrix) <- variables_mlvar_control
  colnames(adj_matrix) <- variables_mlvar_control
  for (i in 1:nrow(participant_data)) {
    predictor <- participant_data$predictor[i]
    criterion <- participant_data$criterion[i]
    relationship <- participant_data$est[i]
    adj_matrix[predictor, criterion] <- relationship
  }
  return(adj_matrix)
}

# Group data by participant and create adjacency matrices
adjacency_matrices_mlvar_control <- results_mlvar_control %>%
  group_by(Participant) %>%
  group_split() %>%
  lapply(create_adjacency_matrix)

# Naming the list elements by Participant IDs 
names(adjacency_matrices_mlvar_control) <- unique(results_mlvar_control$Participant)

# ---------------------------------------------------------------------------- #
# Extracting the person-specific parameters for the 7-node ML-VAR network with fun ----
# ---------------------------------------------------------------------------- #

parameters_mlvar_solution_fun <- mlvarfit_fun$parameters$wilevel.standardized$stdyx.standardized

# ---------------------------------------------------------------------------- #
# Creating the adjacency matrix from ML-VAR parameters for the 7-node network with fun ----
# ---------------------------------------------------------------------------- #

mlvarfit_fun <- readRDS(file = "./02_networks/results/mlvarfit_fun.RDS")

results_mlvar_fun <- parameters_mlvar_solution_fun

# Subset the dataframe to keep only the desired rows
results_mlvar_fun <- results_mlvar_fun %>%
  group_by(cluster) %>%
  slice_head(n = 49) %>%
  ungroup()

#create a column in the dataframe that signifies to which participant do the given set of coefficients belong. 
# Calculate the number of participants
num_participants <- ceiling(nrow(results_mlvar_fun) / 49)

# Create the participant labels
participant_labels <- rep(paste0("Participant ", 1:num_participants), each = 49, length.out = nrow(results_mlvar_fun))

# Add the new column to the dataframe
results_mlvar_fun$Participant <- participant_labels

results_mlvar_fun <- results_mlvar_fun[, c("Participant", "paramHeader", "param", "est", "lower_2.5ci", "upper_2.5ci", "sig", "cluster")]

#create the column with criterion variable names

variables <- c("bad", "energy", "focus", "fun", "interest", "movement", "sad")

repeated_variables <- rep(rep(variables, each = 7), length.out = nrow(results_mlvar_fun))

results_mlvar_fun$criterion <- repeated_variables

#create the column with predictor variable names
repeated_predictor_variables <- rep(variables, length.out = nrow(results_mlvar_fun))

results_mlvar_fun$predictor <- repeated_predictor_variables

#replacing non-significant edges with zeroes

results_mlvar_fun$est <- ifelse((results_mlvar_fun$lower_2.5ci > 0 & results_mlvar_fun$upper_2.5ci > 0) |
                                      (results_mlvar_fun$lower_2.5ci < 0 & results_mlvar_fun$upper_2.5ci < 0), results_mlvar_fun$est, 0)

# ---------------------------------------------------------------------------- #
# Creating the adjacency matrix from ML-VAR parameters for the 7-node network with fun ----
# ---------------------------------------------------------------------------- #

#create an adjacency matrix of auto-regressive and cross-regressive coefficients
#create a directed adjacency matrix for each participant

# Function to create adjacency matrix for a given participant's data
create_adjacency_matrix <- function(participant_data) {
  variables_mlvar_fun <- unique(c(participant_data$criterion, participant_data$predictor))
  adj_matrix <- matrix(NA, nrow = length(variables_mlvar_fun), ncol = length(variables_mlvar_fun))
  rownames(adj_matrix) <- variables_mlvar_fun
  colnames(adj_matrix) <- variables_mlvar_fun
  for (i in 1:nrow(participant_data)) {
    predictor <- participant_data$predictor[i]
    criterion <- participant_data$criterion[i]
    relationship <- participant_data$est[i]
    adj_matrix[predictor, criterion] <- relationship
  }
  return(adj_matrix)
}

# Group data by participant and create adjacency matrices
adjacency_matrices_mlvar_fun <- results_mlvar_fun %>%
  group_by(Participant) %>%
  group_split() %>%
  lapply(create_adjacency_matrix)

# Naming the list elements by Participant IDs 
names(adjacency_matrices_mlvar_fun) <- unique(results_mlvar_fun$Participant)

# ---------------------------------------------------------------------------- #
# TODO: Export temporal network results ----
# ---------------------------------------------------------------------------- #




