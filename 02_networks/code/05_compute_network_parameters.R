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

# # TODO (Determine which ones are needed): Load packages

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

# TODO: Condense code below using custom functions




# ---------------------------------------------------------------------------- #
# Compute expected influence centrality for temporal VAR and ML-VAR networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute outgoing one- and two-step expected influence centralities
# (note that these include both autoregressive and cross-lagged effects)

compute_exp_inf_cent <- function(adj_mats) {
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
  
  return(exp_inf_cent_df)
}

# Run function

exp_inf_cent_var           <- compute_exp_inf_cent(thres_adj_mats_var)
exp_inf_cent_var_control   <- compute_exp_inf_cent(thres_adj_mats_var_control)
exp_inf_cent_var_fun       <- compute_exp_inf_cent(thres_adj_mats_var_fun)

exp_inf_cent_mlvar         <- compute_exp_inf_cent(thres_adj_mats_mlvar)
exp_inf_cent_mlvar_control <- compute_exp_inf_cent(thres_adj_mats_mlvar_control)
exp_inf_cent_mlvar_fun     <- compute_exp_inf_cent(thres_adj_mats_mlvar_fun)





# TODO: Address the following code for 8-node VAR network

all_person_spec_parameters <- expected_influence_wide

#compute mean expected influence of the 2 core depression symptoms (sad and interest)
all_person_spec_parameters <- all_person_spec_parameters %>%
  rowwise() %>%
  mutate(m_ei_s1_sad_int = mean(c(sad_step1, interest_step1), na.rm = TRUE),
         m_ei_s2_sad_int = mean(c(sad_step2, interest_step2), na.rm = TRUE)) %>%
  ungroup()

# TODO: Address the following code for 8-node ML-VAR network

#compute mean expected influence of the 2 core depression symptoms (sad and interest)
all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  rowwise() %>%
  mutate(m_ei_s1_sad_int = mean(c(sad_step1, interest_step1),na.rm = TRUE),
         m_ei_s2_sad_int = mean(c(sad_step2, interest_step2), na.rm = TRUE)) %>%
  ungroup()





# ---------------------------------------------------------------------------- #
# Compute average controllability centrality for temporal VAR and ML-VAR networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute average controllability centrality

compute_avg_cont_cent <- function(adj_mats) {
  # Obtain node names
  
  node_vars <- colnames(adj_mats[[1]])
  
  # Compute average controllability centrality and name result for each node
  
  avg_cont_cent <- lapply(adj_mats, function(x) ave_control_centrality(x))
  
  avg_cont_cent <- lapply(avg_cont_cent, function(x) {
    names(x) <- paste0(node_vars, "_acc")
    
    return(x)
  })
  
  # Transform resulting list into data frame
  
  avg_cont_cent <- as.data.frame(do.call(rbind, avg_cont_cent))
  
  avg_cont_cent$lifepak_id <- rownames(avg_cont_cent)
  
  rownames(avg_cont_cent) <- NULL
  
  avg_cont_cent <- avg_cont_cent %>% relocate(lifepak_id)
  
  return(avg_cont_cent)
}

# Run function

avg_cont_cent_var           <- compute_avg_cont_cent(thres_adj_mats_var)
avg_cont_cent_var_control   <- compute_avg_cont_cent(thres_adj_mats_var_control)
avg_cont_cent_var_fun       <- compute_avg_cont_cent(thres_adj_mats_var_fun)

avg_cont_cent_mlvar         <- compute_avg_cont_cent(thres_adj_mats_mlvar)
avg_cont_cent_mlvar_control <- compute_avg_cont_cent(thres_adj_mats_mlvar_control)
avg_cont_cent_mlvar_fun     <- compute_avg_cont_cent(thres_adj_mats_mlvar_fun)





# TODO: Address the following code for 8-node VAR network

#compute mean average controllability centrality of the 2 core depression symptoms (sad and interest)
results_ave_cont_centrality <- results_ave_cont_centrality %>%
  rowwise() %>%
  mutate(m_acc_sad_int = mean(c(sad_acc, interest_acc), na.rm = TRUE)) %>%
  ungroup()

# TODO: Address the following code for 8-node ML-VAR network

#compute mean average controllability centrality of the 2 core depression symptoms (sad and interest)
results_ave_cont_centrality_mlvar <- results_ave_cont_centrality_mlvar %>%
  rowwise() %>%
  mutate(m_acc_sad_int = mean(c(sad_acc, interest_acc), na.rm = TRUE)) %>%
  ungroup()





# ---------------------------------------------------------------------------- #
# Compute modal controllability centrality for temporal VAR and ML-VAR networks ----
# ---------------------------------------------------------------------------- #

# Define function to compute modal controllability centrality

compute_mod_cont_cent <- function(adj_mats) {
  # Obtain node names
  
  node_vars <- colnames(adj_mats[[1]])
  
  # Compute modal controllability centrality and name result for each node
  
  mod_cont_cent <- lapply(adj_mats, function(x) modal_control_centrality(x))
  
  mod_cont_cent <- lapply(mod_cont_cent, function(x) {
    names(x) <- paste0(node_vars, "_mcc")

    return(x)
  })

  # Transform resulting list into data frame

  mod_cont_cent <- as.data.frame(do.call(rbind, mod_cont_cent))

  mod_cont_cent$lifepak_id <- rownames(mod_cont_cent)

  rownames(mod_cont_cent) <- NULL

  mod_cont_cent <- mod_cont_cent %>% relocate(lifepak_id)
  
  return(mod_cont_cent)
}

# Run function

mod_cont_cent_var           <- compute_mod_cont_cent(thres_adj_mats_var)
mod_cont_cent_var_control   <- compute_mod_cont_cent(thres_adj_mats_var_control)
mod_cont_cent_var_fun       <- compute_mod_cont_cent(thres_adj_mats_var_fun)

mod_cont_cent_mlvar         <- compute_mod_cont_cent(thres_adj_mats_mlvar)
mod_cont_cent_mlvar_control <- compute_mod_cont_cent(thres_adj_mats_mlvar_control)
mod_cont_cent_mlvar_fun     <- compute_mod_cont_cent(thres_adj_mats_mlvar_fun)

# TODO: Address the following code for 8-node VAR network

#compute mean modal controllability centrality of the 2 core depression symptoms (sad and interest)
results_modal_cont_centrality <- results_modal_cont_centrality %>%
  rowwise() %>%
  mutate(m_mcc_sad_int = mean(c(sad_mcc, interest_mcc), na.rm = TRUE)) %>%
  ungroup()

#merge the data frame with expected influence & other parameters with data frames with average control centrality 
#and modal control centrality in a single data frame
all_person_spec_parameters <- all_person_spec_parameters %>%
  left_join(results_ave_cont_centrality, by = "Participant") %>%
  left_join(results_modal_cont_centrality, by = "Participant")

# TODO: Address the following code for 8-node ML-VAR network

#compute mean modal controllability centrality of the 2 core depression symptoms (sad and interest)
results_modal_cont_centrality_mlvar <- results_modal_cont_centrality_mlvar %>%
  rowwise() %>%
  mutate(m_mcc_sad_int = mean(c(sad_mcc, interest_mcc), na.rm = TRUE)) %>%
  ungroup()

#merge the data frame with expected influence & other parameters with data frames with average control centrality and 
#modal control centrality in a single data frame
all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  left_join(results_ave_cont_centrality_mlvar, by = "Participant") %>%
  left_join(results_modal_cont_centrality_mlvar, by = "Participant")





# ---------------------------------------------------------------------------- #
# Computing sums of expected influences ----
# ---------------------------------------------------------------------------- #

# VAR

  # Sum of signed outgoing edges connecting perceived agency ("control") to the 2 core depression symptoms ("sad", "interest") for 8-node VAR network

sum_outgoing <- results %>%
  filter(Predictor == "control" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_cont_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

  # Sum of signed outgoing edges connecting positive activity engagement ("fun") to the 2 core depression symptoms ("sad", "interest")

sum_outgoing_1 <- results %>%
  filter(Predictor == "fun" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_fun_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

  #join these two dataframes with rest of parameters

all_person_spec_parameters <- all_person_spec_parameters %>%
  left_join(sum_outgoing, by = "Participant") %>%
  left_join(sum_outgoing_1, by = "Participant")

# ML-VAR

  # Sum of signed outgoing edges connecting perceived agency ("control") to the 2 core depression symptoms ("sad", "interest") for 8-node ML-VAR network

sum_outgoing_mlvar <- results_mlvar %>%
  filter(Predictor == "control" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_cont_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

  # Sum of signed outgoing edges connecting positive activity engagement ("fun") to the 2 core depression symptoms ("sad", "interest")

sum_outgoing_1_mlvar <- results_mlvar %>%
  filter(Predictor == "fun" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_fun_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

  #join these two dataframes with rest of parameters

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  left_join(sum_outgoing_mlvar, by = "Participant") %>%
  left_join(sum_outgoing_1_mlvar, by = "Participant")





# ---------------------------------------------------------------------------- #
# Computing density ----
# ---------------------------------------------------------------------------- #

# VAR (8-node)

  #inter-node connectivity 

sum_off_diag <- function(mat) {
  sum(mat) - sum(diag(mat))
}

inter_node_connectivity <- sapply(adjacency_matrices, sum_off_diag)

all_person_spec_parameters$`inter-node connectivity` <- inter_node_connectivity

  #intra-node connectivity

sum_diag <- function(mat) {
  sum(diag(mat))
}

all_person_spec_parameters$intra_node_connectivity <- sapply(adjacency_matrices, sum_diag)

  #overall connectivity

all_person_spec_parameters$overall_connectivity <- sapply(adjacency_matrices, function(mat) sum(mat))

# ML-VAR (8-node)

  #inter-node connectivity 

sum_off_diag <- function(mat) {
  sum(mat) - sum(diag(mat))
}

inter_node_connectivity_mlvar <- sapply(adjacency_matrices_mlvar, sum_off_diag)

all_person_spec_parameters_mlvar$`inter-node connectivity` <- inter_node_connectivity_mlvar

  #intra-node connectivity

sum_diag <- function(mat) {
  sum(diag(mat))
}

all_person_spec_parameters_mlvar$intra_node_connectivity <- sapply(adjacency_matrices_mlvar, sum_diag)

  #overall connectivity

all_person_spec_parameters_mlvar$overall_connectivity <- sapply(adjacency_matrices_mlvar, function(mat) sum(mat))






# TODO: Address the following code re expected influence for 7-node VAR network with control

expected_influence_just_control <- expected_influence_control_wide [, c("Participant", "control_step1", "control_step2")]

colnames(expected_influence_just_control) <- c("Participant", "control_step1_7", "control_step2_7")





# TODO: Address the following code re average controllability centrality for 7-node VAR network with control

results_ave_cont_centrality_just_control <- results_ave_cont_centrality_control[, c("Participant", "control_acc")]

colnames(results_ave_cont_centrality_just_control) <- c("Participant", "control_acc_7")





# TODO: Address the following code re modal controllability centrality for 7-node VAR network with control

results_modal_cont_centrality_just_control <- results_modal_cont_centrality_control[, c("Participant", "control_mcc")]

colnames(results_modal_cont_centrality_just_control) <- c("Participant", "control_mcc_7")





# TODO: Address the following code re expected influence for 7-node VAR network with fun

expected_influence_just_fun <- expected_influence_fun_wide [, c("Participant", "fun_step1", "fun_step2")]

colnames (expected_influence_just_fun) <- c("Participant", "fun_step1_7", "fun_step2_7")





# TODO: Address the following code re average controllability centrality for 7-node VAR network with fun

results_ave_cont_centrality_just_fun <- results_ave_cont_centrality_fun [, c("Participant", "fun_acc")]

colnames(results_ave_cont_centrality_just_fun) <- c("Participant", "fun_acc_7")





# TODO: Address the following code re modal controllability centrality for 7-node VAR network with fun

results_modal_cont_centrality_just_fun <- results_modal_cont_centrality_fun[, c("Participant", "fun_mcc")]

colnames(results_modal_cont_centrality_just_fun) <- c("Participant", "fun_mcc_7")





# ---------------------------------------------------------------------------- #
# Merging the entire data frame with person-specific VAR parameters ----
# ---------------------------------------------------------------------------- #

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
# Computing person-specific parameters from the 8-node ML-VAR network ----
# ---------------------------------------------------------------------------- #

# TODO: Address the following code re expected influence for 7-node ML-VAR network with control

expected_influence_just_mlvar_control <- expected_influence_mlvar_control_wide [, c("Participant", "control_step1", "control_step2")]

colnames(expected_influence_just_mlvar_control) <- c("Participant", "control_step1_7", "control_step2_7")




# TODO: Address the following code re average controllability centrality for 7-node ML-VAR network with control

results_ave_cont_centrality_just_mlvar_control <- results_ave_cont_centrality_mlvar_control[, c("Participant", "control_acc")]

colnames(results_ave_cont_centrality_just_mlvar_control) <- c("Participant", "control_acc_7")




# TODO: Address the following code re modal controllability centrality for 7-node ML-VAR network with control

results_modal_cont_centrality_just_mlvar_control <- results_modal_cont_centrality_mlvar_control[, c("Participant", "control_mcc")]

colnames(results_modal_cont_centrality_just_mlvar_control) <- c("Participant", "control_mcc_7")




# TODO: Address the following code re expected influence for 7-node ML-VAR network with fun

expected_influence_just_mlvar_fun <- expected_influence_mlvar_fun_wide [, c("Participant", "fun_step1", "fun_step2")]

colnames(expected_influence_just_mlvar_fun) <- c("Participant", "fun_step1_7", "fun_step2_7")




# TODO: Address the following code re average controllability centrality for 7-node ML-VAR network with fun

results_ave_cont_centrality_just_mlvar_fun <- results_ave_cont_centrality_mlvar_fun [, c("Participant", "fun_acc")]

colnames(results_ave_cont_centrality_just_mlvar_fun) <- c("Participant", "fun_acc_7")





# TODO: Address the following code re modal controllability centrality for 7-node ML-VAR network with fun

results_modal_cont_centrality_just_mlvar_fun <- results_modal_cont_centrality_mlvar_fun[, c("Participant", "fun_mcc")]

colnames(results_modal_cont_centrality_just_mlvar_fun) <- c("Participant", "fun_mcc_7")




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

# TODO: Save parameters as RDS file