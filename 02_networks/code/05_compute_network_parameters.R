# ---------------------------------------------------------------------------- #
# Compute Person-Specific Network Parameters -----
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

# TODO (Determine which ones are needed): Load packages

pkgs <- c("dplyr", "tidyr", "Hmisc", "networktools", "netcontrol")

groundhog.library(pkgs, groundhog_day)





# ---------------------------------------------------------------------------- #
# Import extracted results and thresholded adjacency matrices  ----
# ---------------------------------------------------------------------------- #

extracted_results_path <- "./02_networks/results/extracted/"

load(paste0(extracted_results_path, "results_var.RDS"))
load(paste0(extracted_results_path, "results_var_control.RDS"))
load(paste0(extracted_results_path, "results_var_fun.RDS"))

load(paste0(extracted_results_path, "results_mlvar.RDS"))
load(paste0(extracted_results_path, "results_mlvar_control.RDS"))
load(paste0(extracted_results_path, "results_mlvar_fun.RDS"))

load(paste0(extracted_results_path, "thres_adj_mats_var.RDS"))
load(paste0(extracted_results_path, "thres_adj_mats_var_control.RDS"))
load(paste0(extracted_results_path, "thres_adj_mats_var_fun.RDS"))

load(paste0(extracted_results_path, "thres_adj_mats_mlvar.RDS"))
load(paste0(extracted_results_path, "thres_adj_mats_mlvar_control.RDS"))
load(paste0(extracted_results_path, "thres_adj_mats_mlvar_fun.RDS"))

# ---------------------------------------------------------------------------- #
# Define function used throughout script  ----
# ---------------------------------------------------------------------------- #

# Define function for 7-node networks that restricts a given centrality metric
# to only that for "control" or "fun" (depending on which node was included in
# the network) and labels the centralities as derived from 7-node networks

restrict_to_control_fun_and_label_7 <- function(df, network_type, metric) {
  if (network_type %in% c("7-node with control", "7-node with fun")) {
    if (network_type == "7-node with control") {
      retain_node_vars <- paste0("control_", metric)
    } else if (network_type == "7-node with fun") {
      retain_node_vars <- paste0("fun_",     metric)
    }
    
    df <- df[, c("lifepak_id", retain_node_vars)]
    
    target_cols <- names(df)[names(df) %in% retain_node_vars]
    
    names(df)[names(df) %in% target_cols] <- paste0(target_cols, "_7")
  }
  
  return(df)
}

# ---------------------------------------------------------------------------- #
# Compute expected influence centrality for temporal VAR and ML-VAR networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute outgoing one- and two-step expected influence centralities
# (note that these include both autoregressive and cross-lagged effects)

compute_exp_inf_cent <- function(adj_mats, network_type) {
  # Compute outgoing one- and two-step expected influence centralities
  
  exp_inf_cent <- lapply(adj_mats, function(x) {
    expectedInf(x, step = "both", directed = TRUE)
  })
  
  # Transform resulting nested list into data frame
  
  exp_inf_cent_df <- data.frame()
  
  for (lifepak_id in names(exp_inf_cent)) {
    participant_results <- exp_inf_cent[[lifepak_id]]
    
    for (step in names(participant_results)) {
      step_results <- participant_results[[step]]
      
      temp_df <- data.frame(lifepak_id = lifepak_id,
                            step       = step,
                            t(as.data.frame(step_results)))
      
      exp_inf_cent_df <- rbind(exp_inf_cent_df, temp_df)
    }
  }
  
  rownames(exp_inf_cent_df) <- NULL
  
  # Transform to wide format
  
  possible_node_vars <- c("bad", "control", "energy", "focus", "fun", "interest", "movement", "sad")
  
  node_vars <- names(exp_inf_cent_df)[names(exp_inf_cent_df) %in% possible_node_vars]
  
  exp_inf_cent_df <- exp_inf_cent_df %>%
    pivot_wider(names_from  = step,
                values_from = all_of(node_vars),
                names_sep   = "_")
  
  exp_inf_cent_df <- exp_inf_cent_df[, c("lifepak_id",
                                         names(exp_inf_cent_df)[grep("_step1", names(exp_inf_cent_df))],
                                         names(exp_inf_cent_df)[grep("_step2", names(exp_inf_cent_df))])]
  
  # For 7-node networks, restrict to expected influence centralities of only
  # "control" or "fun" (depending on which node was included in the network)
  # and label the centralities as derived from 7-node networks
  
  metric <- c("step1", "step2")
  
  exp_inf_cent_df <- restrict_to_control_fun_and_label_7(exp_inf_cent_df, network_type, metric)
  
  return(exp_inf_cent_df)
}

# Run function

exp_inf_cent_var           <- compute_exp_inf_cent(thres_adj_mats_var,           "8-node")
exp_inf_cent_var_control   <- compute_exp_inf_cent(thres_adj_mats_var_control,   "7-node with control")
exp_inf_cent_var_fun       <- compute_exp_inf_cent(thres_adj_mats_var_fun,       "7-node with fun")

exp_inf_cent_mlvar         <- compute_exp_inf_cent(thres_adj_mats_mlvar,         "8-node")
exp_inf_cent_mlvar_control <- compute_exp_inf_cent(thres_adj_mats_mlvar_control, "7-node with control")
exp_inf_cent_mlvar_fun     <- compute_exp_inf_cent(thres_adj_mats_mlvar_fun,     "7-node with fun")

# ---------------------------------------------------------------------------- #
# Compute average controllability centrality for temporal VAR and ML-VAR networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute average controllability centrality

compute_avg_cont_cent <- function(adj_mats, network_type) {
  # Obtain node names
  
  node_vars <- colnames(adj_mats[[1]])
  
  # Compute average controllability centrality and name result for each node
  
  avg_cont_cent <- lapply(adj_mats, function(x) ave_control_centrality(x))
  
  avg_cont_cent <- lapply(avg_cont_cent, function(x) {
    names(x) <- paste0(node_vars, "_acc")
    
    return(x)
  })
  
  # Transform resulting list into data frame
  
  avg_cont_cent_df <- as.data.frame(do.call(rbind, avg_cont_cent))
  
  avg_cont_cent_df$lifepak_id <- rownames(avg_cont_cent_df)
  
  rownames(avg_cont_cent_df) <- NULL
  
  avg_cont_cent_df <- avg_cont_cent_df %>% relocate(lifepak_id)
  
  # For 7-node networks, restrict to average controllability centralities of only
  # "control" or "fun" (depending on which node was included in the network)
  # and label the centralities as derived from 7-node networks
  
  metric <- "acc"
  
  avg_cont_cent_df <- restrict_to_control_fun_and_label_7(avg_cont_cent_df, network_type, metric)
  
  return(avg_cont_cent_df)
}

# Run function

avg_cont_cent_var           <- compute_avg_cont_cent(thres_adj_mats_var,           "8-node")
avg_cont_cent_var_control   <- compute_avg_cont_cent(thres_adj_mats_var_control,   "7-node with control")
avg_cont_cent_var_fun       <- compute_avg_cont_cent(thres_adj_mats_var_fun,       "7-node with fun")

avg_cont_cent_mlvar         <- compute_avg_cont_cent(thres_adj_mats_mlvar,         "8-node")
avg_cont_cent_mlvar_control <- compute_avg_cont_cent(thres_adj_mats_mlvar_control, "7-node with control")
avg_cont_cent_mlvar_fun     <- compute_avg_cont_cent(thres_adj_mats_mlvar_fun,     "7-node with fun")

# ---------------------------------------------------------------------------- #
# Compute modal controllability centrality for temporal VAR and ML-VAR networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute modal controllability centrality

compute_mod_cont_cent <- function(adj_mats, network_type) {
  # Obtain node names
  
  node_vars <- colnames(adj_mats[[1]])
  
  # Compute modal controllability centrality and name result for each node
  
  mod_cont_cent <- lapply(adj_mats, function(x) modal_control_centrality(x))
  
  mod_cont_cent <- lapply(mod_cont_cent, function(x) {
    names(x) <- paste0(node_vars, "_mcc")

    return(x)
  })

  # Transform resulting list into data frame

  mod_cont_cent_df <- as.data.frame(do.call(rbind, mod_cont_cent))

  mod_cont_cent_df$lifepak_id <- rownames(mod_cont_cent_df)

  rownames(mod_cont_cent_df) <- NULL

  mod_cont_cent_df <- mod_cont_cent_df %>% relocate(lifepak_id)
  
  # For 7-node networks, restrict to modal controllability centralities of only
  # "control" or "fun" (depending on which node was included in the network)
  # and label the centralities as derived from 7-node networks
  
  metric <- "mcc"
  
  mod_cont_cent_df <- restrict_to_control_fun_and_label_7(mod_cont_cent_df, network_type, metric)
  
  return(mod_cont_cent_df)
}

# Run function

mod_cont_cent_var           <- compute_mod_cont_cent(thres_adj_mats_var,           "8-node")
mod_cont_cent_var_control   <- compute_mod_cont_cent(thres_adj_mats_var_control,   "7-node with control")
mod_cont_cent_var_fun       <- compute_mod_cont_cent(thres_adj_mats_var_fun,       "7-node with fun")

mod_cont_cent_mlvar         <- compute_mod_cont_cent(thres_adj_mats_mlvar,         "8-node")
mod_cont_cent_mlvar_control <- compute_mod_cont_cent(thres_adj_mats_mlvar_control, "7-node with control")
mod_cont_cent_mlvar_fun     <- compute_mod_cont_cent(thres_adj_mats_mlvar_fun,     "7-node with fun")

# ---------------------------------------------------------------------------- #
# Compute mean centralities of two core depression symptoms in 8-node networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute mean of individual centralities for "sad" and "interest"
# in 8-node VAR and ML-VAR networks

compute_m_sad_interest <- function(cent_df, metrics) {
  for (i in 1:length(metrics)) {
    metric <- metrics[i]
    
    sad_name            <- paste0("sad_",            metric)
    interest_name       <- paste0("interest_",       metric)
    m_sad_interest_name <- paste0("m_sad_interest_", metric)
    
    cent_df[, m_sad_interest_name] <- NA
    
    cent_df[, m_sad_interest_name] <- rowMeans(cent_df[, c(sad_name, interest_name)], na.rm = TRUE)
  }
  
  return(cent_df)
}

# Run function for 8-node VAR and ML-VAR networks

exp_inf_cent_var    <- compute_m_sad_interest(exp_inf_cent_var,    c("step1", "step2"))
exp_inf_cent_mlvar  <- compute_m_sad_interest(exp_inf_cent_mlvar,  c("step1", "step2"))

avg_cont_cent_var   <- compute_m_sad_interest(avg_cont_cent_var,   "acc")
avg_cont_cent_mlvar <- compute_m_sad_interest(avg_cont_cent_mlvar, "acc")

mod_cont_cent_var   <- compute_m_sad_interest(mod_cont_cent_var,   "mcc")
mod_cont_cent_mlvar <- compute_m_sad_interest(mod_cont_cent_mlvar, "mcc")

# ---------------------------------------------------------------------------- #
# Compute density in 8-node networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute global expected influence in terms of internode, 
# intranode, and overall connectivity in 8-node VAR and ML-VAR networks

compute_global_exp_inf <- function(adj_mats) {
  # Compute internode, intranode, and overall connectivity
  
  inter_conn <- sapply(adj_mats, function(x) {
    sum(c(x[upper.tri(x)],
          x[lower.tri(x)]))
  })
  
  intra_conn <- sapply(adj_mats, function(x) {
    sum(diag(x))
  })
  
  overall_conn <- sapply(adj_mats, function(x) {
    sum(x)
  })
  
  # Combine into data frame
  
  global_exp_inf_df <- data.frame(lifepak_id   = names(adj_mats),
                                  inter_conn   = inter_conn,
                                  intra_conn   = intra_conn,
                                  overall_conn = overall_conn)
  
  return(global_exp_inf_df)
}

# Run function for 8-node VAR and ML-VAR networks

global_exp_inf_var   <- compute_global_exp_inf(thres_adj_mats_var)
global_exp_inf_mlvar <- compute_global_exp_inf(thres_adj_mats_mlvar)

# ---------------------------------------------------------------------------- #
# Compute sums of selected one-step expected influences in 8-node networks  ----
# ---------------------------------------------------------------------------- #

# Define function to compute (a) sum of signed outgoing edges connecting "control"
# to two core depression symptoms ("sad" and "interest") and (b) same connecting 
# "fun" to such symptoms in 8-node VAR and ML-VAR networks

compute_select_exp_inf <- function(adj_mats) {
  # Compute selected one-step expected influences
  
  sum_control_to_sad_interest <- sapply(adj_mats, function(x) {
    sum(x["control", c("sad", "interest")])
  })
  
  sum_fun_to_sad_interest     <- sapply(adj_mats, function(x) {
    sum(x["fun",     c("sad", "interest")])
  })
  
  # Combine into data frame
  
  select_exp_inf_df <- data.frame(lifepak_id                  = names(adj_mats),
                                  sum_control_to_sad_interest = sum_control_to_sad_interest,
                                  sum_fun_to_sad_interest     = sum_fun_to_sad_interest)
  
  return(select_exp_inf_df)
}

# Run function for 8-node VAR and ML-VAR networks

select_exp_inf_var   <- compute_select_exp_inf(thres_adj_mats_var)
select_exp_inf_mlvar <- compute_select_exp_inf(thres_adj_mats_mlvar)

# ---------------------------------------------------------------------------- #
# Label network parameters and merge into one data frame ----
# ---------------------------------------------------------------------------- #

# Define function to label network parameters and merge into one data frame

merge_net_params <- function(exp_inf_cent,  exp_inf_cent_control,  exp_inf_cent_fun,
                             avg_cont_cent, avg_cont_cent_control, avg_cont_cent_fun,
                             mod_cont_cent, mod_cont_cent_control, mod_cont_cent_fun,
                             global_exp_inf,
                             select_exp_inf) {
  # TODO: Label network parameters
  
  
  
  
  
  # TODO: Append column names and labels for ML-VAR
  
  
  
  
  
  # Merge into one data frame
  
  dfs <- list(exp_inf_cent,  exp_inf_cent_control,  exp_inf_cent_fun,
              avg_cont_cent, avg_cont_cent_control, avg_cont_cent_fun,
              mod_cont_cent, mod_cont_cent_control, mod_cont_cent_fun,
              global_exp_inf,
              select_exp_inf)
  
  merged_df <- Reduce(function(x, y) merge(x, y, by = "lifepak_id", all = TRUE), 
                      dfs)
  
  return(merged_df)
}

# (TODO: Check the merging) Run function for VAR and ML-VAR network parameters

net_params_var   <- merge_net_params(exp_inf_cent_var,    exp_inf_cent_var_control,    exp_inf_cent_var_fun, 
                                     avg_cont_cent_var,   avg_cont_cent_var_control,   avg_cont_cent_var_fun, 
                                     mod_cont_cent_var,   mod_cont_cent_var_control,   mod_cont_cent_var_fun, 
                                     global_exp_inf_var, 
                                     select_exp_inf_var)

net_params_mlvar <- merge_net_params(exp_inf_cent_mlvar,  exp_inf_cent_mlvar_control,  exp_inf_cent_mlvar_fun, 
                                     avg_cont_cent_mlvar, avg_cont_cent_mlvar_control, avg_cont_cent_mlvar_fun, 
                                     mod_cont_cent_mlvar, mod_cont_cent_mlvar_control, mod_cont_cent_mlvar_fun, 
                                     global_exp_inf_mlvar, 
                                     select_exp_inf_mlvar)





# ---------------------------------------------------------------------------- #
# Merging the entire data frame with person-specific VAR parameters ----
# ---------------------------------------------------------------------------- #

#merge the data frame with expected influence & other parameters with data frames with average control centrality 
#and modal control centrality in a single data frame
all_person_spec_parameters <- all_person_spec_parameters %>%
  left_join(results_ave_cont_centrality, by = "Participant") %>%
  left_join(results_modal_cont_centrality, by = "Participant")

#merge the data frame with expected influence & other parameters with data frames with average control centrality and 
#modal control centrality in a single data frame
all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  left_join(results_ave_cont_centrality_mlvar, by = "Participant") %>%
  left_join(results_modal_cont_centrality_mlvar, by = "Participant")

all_person_spec_parameters <- all_person_spec_parameters %>%
  left_join(expected_influence_just_fun, by = "Participant") %>%
  left_join(expected_influence_just_control, by = "Participant") %>%
  left_join(results_ave_cont_centrality_just_control, by = "Participant") %>%
  left_join(results_modal_cont_centrality_just_control, by = "Participant") %>%
  left_join(results_ave_cont_centrality_just_fun, by = "Participant") %>%
  left_join(results_modal_cont_centrality_just_fun, by = "Participant")

# ---------------------------------------------------------------------------- #
# Labeling columns for person-specific VAR parameters ----
# ---------------------------------------------------------------------------- #

all_person_spec_parameters <- all_person_spec_parameters %>%
  mutate_at(vars(bad_step1, control_step1, energy_step1, focus_step1, fun_step1, interest_step1, movement_step1, sad_step1), ~ {
    label(.) <- "One-step expected influence for a given node"
    .
  })

all_person_spec_parameters <- all_person_spec_parameters %>%
  mutate_at(vars(bad_step2, control_step2, energy_step2, focus_step2, fun_step2, interest_step2, movement_step2, sad_step2), ~ {
    label(.) <- "Two-step expected influence for a given node"
    .
  })

label(all_person_spec_parameters$m_ei_s1_sad_int)    <- "Mean one-step expected influence of sad and interest"
label(all_person_spec_parameters$m_ei_s2_sad_int)    <- "Mean two-step expected influence of sad and interest"
label(all_person_spec_parameters$sum_cont_to_sadint) <- "Sum of outgoing edges connecting control to sad and interest"
label(all_person_spec_parameters$sum_fun_to_sadint)  <- "Sum of outgoing edges connecting fun to sad and interest"

all_person_spec_parameters <- all_person_spec_parameters %>%
  mutate_at(vars(bad_acc, control_acc, energy_acc, focus_acc, 
                 fun_acc, interest_acc, movement_acc, sad_acc), ~ {
                   label(.) <- "Average controllability centrality for a given node"
                   .
                 })

label(all_person_spec_parameters$m_acc_sad_int) <- "Mean average controllability centrality of sad and interest"

all_person_spec_parameters <- all_person_spec_parameters %>%
  mutate_at(vars(bad_mcc, control_mcc, energy_mcc, focus_mcc, 
                 fun_mcc, interest_mcc, movement_mcc, sad_mcc), ~ {
                   label(.) <- "Modal controllability centrality for a given node"
                   .
                 })

label(all_person_spec_parameters$m_mcc_sad_int)   <- "Mean modal controllability centrality of sad and interest"
label(all_person_spec_parameters$fun_step1_7)     <- "One-step expected influence of fun in 7-node network"
label(all_person_spec_parameters$fun_step2_7)     <- "Two-step expected influence of fun in 7-node network"
label(all_person_spec_parameters$control_step1_7) <- "One-step expected influence of control in 7-node network"
label(all_person_spec_parameters$control_step2_7) <- "Two-step expected influence of control in 7-node network"
label(all_person_spec_parameters$control_acc_7)   <- "Average controllability centrality of control in a 7-node network"
label(all_person_spec_parameters$fun_acc_7)       <- "Average controllability centrality of fun in a 7-node network"
label(all_person_spec_parameters$control_mcc_7)   <- "Modal controllability centrality of control in a 7-node network"
label(all_person_spec_parameters$fun_mcc_7)       <- "Modal controllability centrality of fun in a 7-node network"

# ---------------------------------------------------------------------------- #
# Merging the entire data frame with person-specific parameters from ML-VAR ----
# ---------------------------------------------------------------------------- #

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  left_join(expected_influence_just_mlvar_fun, by = "Participant") %>%
  left_join(expected_influence_just_mlvar_control, by = "Participant") %>%
  left_join(results_ave_cont_centrality_just_mlvar_control, by = "Participant") %>%
  left_join(results_modal_cont_centrality_just_mlvar_control, by = "Participant") %>%
  left_join(results_ave_cont_centrality_just_mlvar_fun, by = "Participant") %>%
  left_join(results_modal_cont_centrality_just_mlvar_fun, by = "Participant")

# ---------------------------------------------------------------------------- #
# Labeling the columns for person-specific parameters from ML-VAR ----
# ---------------------------------------------------------------------------- #

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  mutate_at(vars(bad_step1, control_step1, energy_step1, focus_step1, fun_step1, interest_step1, movement_step1, sad_step1), ~ {
    label(.) <- "ML_VAR One-step expected influence for a given node"
    .
  })

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  mutate_at(vars(bad_step2, control_step2, energy_step2, focus_step2, fun_step2, interest_step2, movement_step2, sad_step2), ~ {
    label(.) <- "ML_VAR Two-step expected influence for a given node"
    .
  })

label(all_person_spec_parameters_mlvar$m_ei_s1_sad_int)    <- "ML_VAR Mean one-step expected influence of sad and interest"
label (all_person_spec_parameters_mlvar$m_ei_s2_sad_int)   <- "ML_VAR Mean two-step expected influence of sad and interest"
label(all_person_spec_parameters_mlvar$sum_cont_to_sadint) <- "ML_VAR Sum of outgoing edges connecting control to sad and interest"
label(all_person_spec_parameters_mlvar$sum_fun_to_sadint)  <- "ML_VAR Sum of outgoing edges connecting fun to sad and interest"

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  mutate_at(vars(bad_acc, control_acc, energy_acc, focus_acc, 
                 fun_acc, interest_acc, movement_acc, sad_acc), ~ {
                   label(.) <- "ML_VAR Average controllability centrality for a given node"
                   .
                 })

label(all_person_spec_parameters_mlvar$m_acc_sad_int) <- "ML_VAR Mean average controllability centrality of sad and interest"

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  mutate_at(vars(bad_mcc, control_mcc, energy_mcc, focus_mcc, 
                 fun_mcc, interest_mcc, movement_mcc, sad_mcc), ~ {
                   label(.) <- "ML_VAR Modal controllability centrality for a given node"
                   .
                 })

label(all_person_spec_parameters_mlvar$m_mcc_sad_int)   <- "ML_VAR Mean modal controllability centrality of sad and interest"
label(all_person_spec_parameters_mlvar$fun_step1_7)     <- "ML_VAR One-step expected influence of fun in 7-node network"
label(all_person_spec_parameters_mlvar$fun_step2_7)     <- "ML_VAR Two-step expected influence of fun in 7-node network"
label(all_person_spec_parameters_mlvar$control_step1_7) <- "ML_VAR One-step expected influence of control in 7-node network"
label(all_person_spec_parameters_mlvar$control_step2_7) <- "ML_VAR Two-step expected influence of control in 7-node network"
label(all_person_spec_parameters_mlvar$control_acc_7)   <- "ML_VAR Average controllability centrality of control in a 7-node network"
label(all_person_spec_parameters_mlvar$fun_acc_7)       <- "ML_VAR Average controllability centrality of fun in a 7-node network"
label(all_person_spec_parameters_mlvar$control_mcc_7)   <- "ML_VAR Modal controllability centrality of control in a 7-node network"
label(all_person_spec_parameters_mlvar$fun_mcc_7)       <- "ML_VAR Modal controllability centrality of fun in a 7-node network"

# ---------------------------------------------------------------------------- #
# Merging the data frames with person-specific parameters from VAR and person-specific parameters extracted from ML-VAR ----
# ---------------------------------------------------------------------------- #

#changing the names of all columns in the ML-VAR data frame so that they end with _mlvar
colnames(all_person_spec_parameters_mlvar) <- paste0(colnames(all_person_spec_parameters_mlvar), "_mlvar")

colnames(all_person_spec_parameters_mlvar)[colnames(all_person_spec_parameters_mlvar) == "Participant_mlvar"] <- "Participant"

#merging the two data frames with person-specific parameters

all_person_spec_parameters_var_mlvar <- all_person_spec_parameters %>%
  left_join(all_person_spec_parameters_mlvar, by = "Participant")

# ---------------------------------------------------------------------------- #
# Export network parameters ----
# ---------------------------------------------------------------------------- #

# TODO: Save parameters as RDS file




