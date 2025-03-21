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

# Load packages

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

# Define function to label network parameters in given data frame

label_net_params <- function(df, target_cols, labels) {
  for (i in 1:length(target_cols)) {
    target_col <- target_cols[i]
    
    if (length(labels) == 1) {
      label <- labels
    } else if (length(labels) == length(target_cols)) {
      label <- labels[i]
    }
    
    label(df[[target_col]]) <- label
  }
  
  return(df)
}

# Define function to label network parameters and merge into one data frame

label_merge_net_params <- function(exp_inf_cent,  exp_inf_cent_control,  exp_inf_cent_fun,
                                   avg_cont_cent, avg_cont_cent_control, avg_cont_cent_fun,
                                   mod_cont_cent, mod_cont_cent_control, mod_cont_cent_fun,
                                   global_exp_inf, select_exp_inf, model_type) {
  # Label network parameters
  
  node_vars <- c("bad", "control", "energy", "focus", "fun", "interest", "movement", "sad")
  
  exp_inf_cent  <- label_net_params(exp_inf_cent,
                                    paste0(node_vars, "_step1"),
                                    "one-step expected influence centrality for given node in 8-node network")
  exp_inf_cent  <- label_net_params(exp_inf_cent,
                                    paste0(node_vars, "_step2"),
                                    "two-step expected influence centrality for given node in 8-node network")
  
  avg_cont_cent <- label_net_params(avg_cont_cent,
                                    paste0(node_vars, "_acc"),
                                    "average controllability centrality for given node in 8-node network")
  
  mod_cont_cent <- label_net_params(mod_cont_cent,
                                    paste0(node_vars, "_mcc"),
                                    "model controllability centrality for given node in 8-node network")
  
  label(exp_inf_cent$m_sad_interest_step1)          <- "mean of one-step expected influence centralities for sad and interest in 8-node network"
  label(exp_inf_cent$m_sad_interest_step2)          <- "mean of two-step expected influence centralities for sad and interest in 8-node network"
  
  label(avg_cont_cent$m_sad_interest_acc)           <- "mean of average controllability centralities for sad and interest in 8-node network"
  
  label(mod_cont_cent$m_sad_interest_mcc)           <- "mean of modal controllability centralities for sad and interest in 8-node network"
  
  label(exp_inf_cent_control$control_step1_7)       <- "one-step expected influence centrality for control in 7-node network"
  label(exp_inf_cent_control$control_step2_7)       <- "two-step expected influence centrality for control in 7-node network"
    
  label(exp_inf_cent_fun$fun_step1_7)               <- "one-step expected influence centrality for fun in 7-node network"
  label(exp_inf_cent_fun$fun_step2_7)               <- "two-step expected influence centrality for fun in 7-node network"
  
  label(avg_cont_cent_control$control_acc_7)        <- "average controllability centrality for control in 7-node network"

  label(avg_cont_cent_fun$fun_acc_7)                <- "average controllability centrality for fun in 7-node network"
  
  label(mod_cont_cent_control$control_mcc_7)        <- "modal controllability centrality for control in 7-node network"
  
  label(mod_cont_cent_fun$fun_mcc_7)                <- "modal controllability centrality for fun in 7-node network"
  
  label(global_exp_inf$inter_conn)                  <- "global internode expected influence in 8-node network"
  label(global_exp_inf$intra_conn)                  <- "global intranode expected influence in 8-node network"
  label(global_exp_inf$overall_conn)                <- "global overall expected influence in 8-node network"
  
  label(select_exp_inf$sum_control_to_sad_interest) <- "sum of signed edges from control to sad and interest in 8-node network"
  label(select_exp_inf$sum_fun_to_sad_interest)     <- "sum of signed edges from fun to sad and interest in 8-node network"
  
  # Merge into one data frame
  
  dfs <- list(exp_inf_cent,  exp_inf_cent_control,  exp_inf_cent_fun,
              avg_cont_cent, avg_cont_cent_control, avg_cont_cent_fun,
              mod_cont_cent, mod_cont_cent_control, mod_cont_cent_fun,
              global_exp_inf,
              select_exp_inf)
  
  merged_df <- Reduce(function(x, y) merge(x, y, by = "lifepak_id", all = TRUE), 
                      dfs)
  
  merged_df$lifepak_id <- as.integer(merged_df$lifepak_id)
  
  merged_df <- merged_df[order(merged_df$lifepak_id), ]
  
  row.names(merged_df) <- NULL
  
  # Append column names and labels based on model type
  
  target_colnames <- names(merged_df)[names(merged_df) != "lifepak_id"]
  target_labels   <- label(merged_df[target_colnames])
  
  if (model_type == "VAR") {
    suffix <- "_var"
  } else if (model_type == "ML-VAR") {
    suffix <- "_mlvar"
  }
  
  new_colnames <- paste0(target_colnames, "_", suffix)
  new_labels   <- paste(model_type, target_labels)
  
  names(merged_df)[names(merged_df) %in% target_colnames] <- new_colnames
  
  merged_df <- label_net_params(merged_df, new_colnames, new_labels)
  
  return(merged_df)
}

# Run function for VAR and ML-VAR network parameters

net_params_var   <- label_merge_net_params(exp_inf_cent_var,     exp_inf_cent_var_control,    exp_inf_cent_var_fun, 
                                           avg_cont_cent_var,    avg_cont_cent_var_control,   avg_cont_cent_var_fun, 
                                           mod_cont_cent_var,    mod_cont_cent_var_control,   mod_cont_cent_var_fun, 
                                           global_exp_inf_var,   select_exp_inf_var, "VAR")

net_params_mlvar <- label_merge_net_params(exp_inf_cent_mlvar,   exp_inf_cent_mlvar_control,  exp_inf_cent_mlvar_fun, 
                                           avg_cont_cent_mlvar,  avg_cont_cent_mlvar_control, avg_cont_cent_mlvar_fun, 
                                           mod_cont_cent_mlvar,  mod_cont_cent_mlvar_control, mod_cont_cent_mlvar_fun, 
                                           global_exp_inf_mlvar, select_exp_inf_mlvar, "ML-VAR")

# ---------------------------------------------------------------------------- #
# Merge VAR and ML-VAR network parameters into one data frame ----
# ---------------------------------------------------------------------------- #

net_params_var_mlvar <- merge(net_params_var, net_params_mlvar, by = "lifepak_id", all = TRUE)

# ---------------------------------------------------------------------------- #
# Exclude mean centralities and overall connectivity ----
# ---------------------------------------------------------------------------- #

# Given that mean centralities are perfectly collinear with individual centralities
# and that overall connectivity is perfectly collinear with internode and intranode
# connectivities, exclude mean centralities and overall connectivity

exclude_params_stem <- c("m_sad_interest_step1", "m_sad_interest_step2", 
                         "m_sad_interest_acc", "m_sad_interest_mcc", "overall_conn")

exclude_params <- c(paste0(exclude_params_stem, "__var"),
                    paste0(exclude_params_stem, "__mlvar"))

net_params_var_mlvar <- net_params_var_mlvar[, !(names(net_params_var_mlvar) %in% exclude_params)]

# ---------------------------------------------------------------------------- #
# Export network parameters ----
# ---------------------------------------------------------------------------- #

net_params_path <- "./02_networks/results/net_params/"

dir.create(net_params_path)

save(net_params_var_mlvar, file = paste0(net_params_path, "net_params_var_mlvar.RDS"))