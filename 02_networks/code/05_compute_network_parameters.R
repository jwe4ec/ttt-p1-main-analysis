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
# Computing expected influence for the 8-node VAR network ----
# ---------------------------------------------------------------------------- #

#Computing expected influence for participant 1
expectedInf(adjacency_matrices[["Participant 1"]], step = c("both", 1, 2), directed = TRUE)

#Computing expected influence for all participants
apply_expectedInf <- function(matrix) {
  expectedInf(matrix, step = c("both", 1, 2), directed = TRUE)
}

results_exp_inf <- lapply(adjacency_matrices, apply_expectedInf)

#transform this list into a data frame

expected_influence <- data.frame()

for (participant in names(results_exp_inf)) {
  participant_data <- results_exp_inf[[participant]]
  
  
  for (step in names(participant_data)) {
    step_data <- participant_data[[step]]
    
    temp_df <- data.frame(
      Participant = participant,
      Step = step,
      t(as.data.frame(step_data))
    )
    
    expected_influence <- rbind(expected_influence, temp_df)
  }
}

# Reset row names
rownames(expected_influence) <- NULL


expected_influence_wide <- expected_influence %>%
  pivot_wider(
    names_from = Step, 
    values_from = c(bad, control, energy, focus, fun, interest, movement, sad),
    names_sep = "_"
  ) 

expected_influence_wide <- expected_influence_wide[, c("Participant", 
                                                       "bad_step1", "control_step1", "energy_step1", "focus_step1", "fun_step1", "interest_step1", "movement_step1", "sad_step1",
                                                       "bad_step2", "control_step2", "energy_step2", "focus_step2", "fun_step2", "interest_step2", "movement_step2", "sad_step2")]

all_person_spec_parameters <- expected_influence_wide

#compute mean expected influence of the 2 core depression symptoms (sad and interest)
all_person_spec_parameters <- all_person_spec_parameters %>%
  rowwise() %>%
  mutate(m_ei_s1_sad_int = mean(c(sad_step1, interest_step1), na.rm = TRUE),
         m_ei_s2_sad_int = mean(c(sad_step2, interest_step2), na.rm = TRUE)) %>%
  ungroup()

# ---------------------------------------------------------------------------- #
# Computing sum of signed outgoing edges connecting perceived agency ("control") to the 2 core depression symptoms ("sad", "interest") for the 8-node VAR network ----
# ---------------------------------------------------------------------------- #

sum_outgoing <- results %>%
  filter(Predictor == "control" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_cont_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

#compute sum of signed outgoing edges connecting positive activity engagement ("fun") to the 2 core depression symptoms ("sad", "interest")

sum_outgoing_1 <- results %>%
  filter(Predictor == "fun" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_fun_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

#join these two dataframes with the sum of signed outgoing edges with the data frame with the rest of parameters

all_person_spec_parameters <- all_person_spec_parameters %>%
  left_join(sum_outgoing, by = "Participant") %>%
  left_join(sum_outgoing_1, by = "Participant")

# ---------------------------------------------------------------------------- #
# Computing density (inter-node, intra-node and overall connectivity - 3 parameters) for the 8-node VAR network ----
# ---------------------------------------------------------------------------- #

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

# ---------------------------------------------------------------------------- #
# Computing controllability centrality for the 8-node VAR network ----
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Computing average controllability centrality for the 8-node VAR network ----
# ---------------------------------------------------------------------------- #

#Computing average control centrality for all participants and transforming the resulting list into a dataframe

results_ave_cont_centrality <- lapply(adjacency_matrices, function(matrix) ave_control_centrality(matrix))

results_ave_cont_centrality <- do.call(rbind, results_ave_cont_centrality)

results_ave_cont_centrality <- as.data.frame(results_ave_cont_centrality)

results_ave_cont_centrality$Participant <- rownames(results_ave_cont_centrality)

rownames(results_ave_cont_centrality) <- NULL

results_ave_cont_centrality <- results_ave_cont_centrality %>% relocate(Participant, .before = V1)

colnames(results_ave_cont_centrality) <- c("Participant", "bad_acc", "control_acc", "energy_acc", "focus_acc", 
                                           "fun_acc", "interest_acc", "movement_acc", "sad_acc")

#compute mean average controllability centrality of the 2 core depression symptoms (sad and interest)
results_ave_cont_centrality <- results_ave_cont_centrality %>%
  rowwise() %>%
  mutate(m_acc_sad_int = mean(c(sad_acc, interest_acc), na.rm = TRUE)) %>%
  ungroup()

# ---------------------------------------------------------------------------- #
# Computing modal controllability centrality for the 8-node VAR network ----
# ---------------------------------------------------------------------------- #

results_modal_cont_centrality <- lapply(adjacency_matrices, function(matrix) modal_control_centrality(matrix))

results_modal_cont_centrality <- do.call(rbind, results_modal_cont_centrality)

results_modal_cont_centrality <- as.data.frame(results_modal_cont_centrality)

results_modal_cont_centrality$Participant <- rownames(results_modal_cont_centrality)

rownames(results_modal_cont_centrality) <- NULL

results_modal_cont_centrality <- results_modal_cont_centrality %>% relocate(Participant, .before = V1)

colnames(results_modal_cont_centrality) <- c("Participant", "bad_mcc", "control_mcc", "energy_mcc", "focus_mcc", 
                                             "fun_mcc", "interest_mcc", "movement_mcc", "sad_mcc")

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

# ---------------------------------------------------------------------------- #
# Computing expected influence for all participants - for the 7-node VAR network with control ----
# ---------------------------------------------------------------------------- #

results_exp_inf_control <- lapply(adjacency_matrices_control, apply_expectedInf)

#transform this list into a data frame

expected_influence_control <- data.frame()

for (participant in names(results_exp_inf_control)) {
  participant_data <- results_exp_inf_control[[participant]]
  
  for (step in names(participant_data)) {
    step_data <- participant_data[[step]]
    
    temp_df <- data.frame(
      Participant = participant,
      Step = step,
      t(as.data.frame(step_data))
    )
    
    expected_influence_control <- rbind(expected_influence_control, temp_df)
  }
}

rownames(expected_influence_control) <- NULL

expected_influence_control_wide <- expected_influence_control %>%
  pivot_wider(
    names_from = Step, 
    values_from = c(bad, control, energy, focus, interest, movement, sad),
    names_sep = "_"
  ) 

expected_influence_control_wide <- expected_influence_control_wide [, c("Participant", 
                                                                        "bad_step1", "control_step1", "energy_step1", "focus_step1", "interest_step1", "movement_step1", "sad_step1",
                                                                        "bad_step2", "control_step2", "energy_step2", "focus_step2", "interest_step2", "movement_step2", "sad_step2")]

expected_influence_just_control <- expected_influence_control_wide [, c("Participant", "control_step1", "control_step2")]

colnames(expected_influence_just_control) <- c("Participant", "control_step1_7", "control_step2_7")

# ---------------------------------------------------------------------------- #
# Computing average controllability centrality for the 7-node VAR network with control ----
# ---------------------------------------------------------------------------- #

#ave_control_centrality function

#computing average control centrality for all participants and transforming the resulting list into a dataframe

results_ave_cont_centrality_control <- lapply(adjacency_matrices_control, function(matrix) ave_control_centrality(matrix))

results_ave_cont_centrality_control <- do.call(rbind, results_ave_cont_centrality_control)

results_ave_cont_centrality_control <- as.data.frame(results_ave_cont_centrality_control)

results_ave_cont_centrality_control$Participant <- rownames(results_ave_cont_centrality_control)

rownames(results_ave_cont_centrality_control) <- NULL

results_ave_cont_centrality_control <- results_ave_cont_centrality_control %>% relocate(Participant, .before = V1)

colnames(results_ave_cont_centrality_control) <- c("Participant", "bad_acc", "control_acc", "energy_acc", "focus_acc", 
                                                   "interest_acc", "movement_acc", "sad_acc")

results_ave_cont_centrality_just_control <- results_ave_cont_centrality_control [, c("Participant", "control_acc")]

colnames(results_ave_cont_centrality_just_control) <- c("Participant", "control_acc_7")

# ---------------------------------------------------------------------------- #
# Computing modal controllability centrality for the 7-node VAR network with control ----
# ---------------------------------------------------------------------------- #

results_modal_cont_centrality_control <- lapply(adjacency_matrices_control, function(matrix) modal_control_centrality(matrix))

results_modal_cont_centrality_control <- do.call(rbind, results_modal_cont_centrality_control)

results_modal_cont_centrality_control <- as.data.frame(results_modal_cont_centrality_control)

results_modal_cont_centrality_control$Participant <- rownames(results_modal_cont_centrality_control)

rownames(results_modal_cont_centrality_control) <- NULL

results_modal_cont_centrality_control <- results_modal_cont_centrality_control %>% relocate(Participant, .before = V1)

colnames(results_modal_cont_centrality_control) <- c("Participant", "bad_mcc", "control_mcc", "energy_mcc", "focus_mcc", 
                                                     "interest_mcc", "movement_mcc", "sad_mcc")

results_modal_cont_centrality_just_control <- results_modal_cont_centrality_control[, c("Participant", "control_mcc")]

colnames(results_modal_cont_centrality_just_control) <- c("Participant", "control_mcc_7")

# ---------------------------------------------------------------------------- #
# Computing the parameters in a 7-node VAR network with fun ----
# ---------------------------------------------------------------------------- #

#computing expected influence for all participants

results_exp_inf_fun <- lapply(adjacency_matrices_fun, apply_expectedInf)

#transform this list into a data frame

expected_influence_fun <- data.frame()

for (participant in names(results_exp_inf_fun)) {
  participant_data <- results_exp_inf_fun[[participant]]
  
  for (step in names(participant_data)) {
    step_data <- participant_data[[step]]
    
    temp_df <- data.frame(
      Participant = participant,
      Step = step,
      t(as.data.frame(step_data))
    )
    
    expected_influence_fun <- rbind(expected_influence_fun, temp_df)
  }
}

rownames(expected_influence_fun) <- NULL

expected_influence_fun_wide <- expected_influence_fun %>%
  pivot_wider(
    names_from = Step, 
    values_from = c(bad, fun, energy, focus, interest, movement, sad),
    names_sep = "_"
  ) 

expected_influence_fun_wide <- expected_influence_fun_wide [, c("Participant", 
                                                                "bad_step1", "energy_step1", "focus_step1", "fun_step1", "interest_step1", "movement_step1", "sad_step1",
                                                                "bad_step2", "energy_step2", "focus_step2", "fun_step2", "interest_step2", "movement_step2", "sad_step2")]

expected_influence_just_fun <- expected_influence_fun_wide [, c("Participant", "fun_step1", "fun_step2")]

colnames (expected_influence_just_fun) <- c("Participant", "fun_step1_7", "fun_step2_7")

# ---------------------------------------------------------------------------- #
# Computing average contolability centrality for the 7-node VAR network with fun ----
# ---------------------------------------------------------------------------- #

#ave_fun_centrality function

#computing average control centrality for all participants and transforming the resulting list into a dataframe

results_ave_cont_centrality_fun <- lapply(adjacency_matrices_fun, function(matrix) ave_control_centrality(matrix))

results_ave_cont_centrality_fun <- do.call(rbind, results_ave_cont_centrality_fun)

results_ave_cont_centrality_fun <- as.data.frame(results_ave_cont_centrality_fun)

results_ave_cont_centrality_fun$Participant <- rownames(results_ave_cont_centrality_fun)

rownames(results_ave_cont_centrality_fun) <- NULL

results_ave_cont_centrality_fun <- results_ave_cont_centrality_fun %>% relocate(Participant, .before = V1)

colnames(results_ave_cont_centrality_fun) <- c("Participant", "bad_acc", "energy_acc", "focus_acc", 
                                               "fun_acc", "interest_acc", "movement_acc", "sad_acc")

results_ave_cont_centrality_just_fun <- results_ave_cont_centrality_fun [, c("Participant", "fun_acc")]

colnames(results_ave_cont_centrality_just_fun) <- c("Participant", "fun_acc_7")

# ---------------------------------------------------------------------------- #
# Computing modal controlability centrality for the 7-node VAR network with fun ----
# ---------------------------------------------------------------------------- #

results_modal_cont_centrality_fun <- lapply(adjacency_matrices_fun, function(matrix) modal_control_centrality(matrix))

results_modal_cont_centrality_fun <- do.call(rbind, results_modal_cont_centrality_fun)

results_modal_cont_centrality_fun <- as.data.frame(results_modal_cont_centrality_fun)

results_modal_cont_centrality_fun$Participant <- rownames(results_modal_cont_centrality_fun)

rownames(results_modal_cont_centrality_fun) <- NULL

results_modal_cont_centrality_fun <- results_modal_cont_centrality_fun %>% relocate(Participant, .before = V1)

colnames(results_modal_cont_centrality_fun) <- c("Participant", "bad_mcc", "energy_mcc", "focus_mcc", 
                                                 "fun_mcc", "interest_mcc", "movement_mcc", "sad_mcc")

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

#Computing expected influence for all participants
apply_expectedInf <- function(matrix) {
  expectedInf(matrix, step = c("both", 1, 2), directed = TRUE)
}

#Computing expected influence for all participants

results_exp_inf_mlvar <- lapply(adjacency_matrices_mlvar, apply_expectedInf)

#transform this list into a data frame

expected_influence_mlvar <- data.frame()

for (participant in names(results_exp_inf_mlvar)) {
  participant_data <- results_exp_inf_mlvar[[participant]]
  
  for (step in names(participant_data)) {
    step_data <- participant_data[[step]]
    
    temp_df <- data.frame(
      Participant = participant,
      Step = step,
      t(as.data.frame(step_data))
    )
    
    expected_influence_mlvar <- rbind(expected_influence_mlvar, temp_df)
  }
}

rownames(expected_influence_mlvar) <- NULL

expected_influence_wide_mlvar <- expected_influence_mlvar %>%
  pivot_wider(
    names_from = Step, 
    values_from = c(bad, control, energy, focus, fun, interest, movement, sad),
    names_sep = "_"
  ) 

expected_influence_wide_mlvar <- expected_influence_wide_mlvar [, c("Participant", 
                                                                    "bad_step1", "control_step1", "energy_step1", "focus_step1", "fun_step1", "interest_step1", "movement_step1", "sad_step1", 
                                                                    "bad_step2", "control_step2", "energy_step2", "focus_step2", "fun_step2", "interest_step2", "movement_step2", "sad_step2")]

all_person_spec_parameters_mlvar <- expected_influence_wide_mlvar

#compute mean expected influence of the 2 core depression symptoms (sad and interest)
all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  rowwise() %>%
  mutate(m_ei_s1_sad_int = mean(c(sad_step1, interest_step1),na.rm = TRUE),
         m_ei_s2_sad_int = mean(c(sad_step2, interest_step2), na.rm = TRUE)) %>%
  ungroup()

# ---------------------------------------------------------------------------- #
# Computing sum of signed outgoing edges connecting perceived agency ("control") to the 2 core depression symptoms ("sad", "interest") for the 8-node ML-VAR network ----
# ---------------------------------------------------------------------------- #

sum_outgoing_mlvar <- results_mlvar %>%
  filter(Predictor == "control" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_cont_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

#compute sum of signed outgoing edges connecting positive activity engagement ("fun") to the 2 core depression symptoms ("sad", "interest")

sum_outgoing_1_mlvar <- results_mlvar %>%
  filter(Predictor == "fun" & (Criterion == "interest" | Criterion == "sad")) %>% 
  group_by(Participant) %>%
  summarise(sum_fun_to_sadint = sum(est, na.rm = TRUE)) %>%
  mutate(Participant_number = as.numeric(gsub(".*[^0-9]", "", Participant))) %>% #have to sort the "Participant" column according to the numeric suffix
  arrange(Participant_number) %>%
  select(-Participant_number)

#join these two dataframes with the sum of signed outgoind edges with the data frame with the rest of parameters

all_person_spec_parameters_mlvar <- all_person_spec_parameters_mlvar %>%
  left_join(sum_outgoing_mlvar, by = "Participant") %>%
  left_join(sum_outgoing_1_mlvar, by = "Participant")

# ---------------------------------------------------------------------------- #
# Computing ML-VAR density (inter-node, intra-node and overall connectivity - 3 parameters) for the 8-node ML-VAR network ----
# ---------------------------------------------------------------------------- #

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

# ---------------------------------------------------------------------------- #
# Computing ML-VAR controllability centrality for the 8-node ML-VAR network ----
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Computing average control centrality for the 8-node ML-VAR network ----
# ---------------------------------------------------------------------------- #

results_ave_cont_centrality_mlvar <- lapply(adjacency_matrices_mlvar, function(matrix) ave_control_centrality(matrix))

results_ave_cont_centrality_mlvar <- do.call(rbind, results_ave_cont_centrality_mlvar)

results_ave_cont_centrality_mlvar <- as.data.frame(results_ave_cont_centrality_mlvar)

results_ave_cont_centrality_mlvar$Participant <- rownames(results_ave_cont_centrality_mlvar)

rownames(results_ave_cont_centrality_mlvar) <- NULL

results_ave_cont_centrality_mlvar <- results_ave_cont_centrality_mlvar %>% relocate(Participant, .before = V1)

colnames(results_ave_cont_centrality_mlvar) <- c("Participant", "bad_acc", "control_acc", "energy_acc", "focus_acc", 
                                                  "fun_acc", "interest_acc", "movement_acc", "sad_acc")

#compute mean average controllability centrality of the 2 core depression symptoms (sad and interest)
results_ave_cont_centrality_mlvar <- results_ave_cont_centrality_mlvar %>%
  rowwise() %>%
  mutate(m_acc_sad_int = mean(c(sad_acc, interest_acc), na.rm = TRUE)) %>%
  ungroup()

# ---------------------------------------------------------------------------- #
# Computing ML-VAR modal controllability centrality for the 8-node ML-VAR network ----
# ---------------------------------------------------------------------------- #

results_modal_cont_centrality_mlvar <- lapply(adjacency_matrices_mlvar, function(matrix) modal_control_centrality(matrix))

results_modal_cont_centrality_mlvar <- do.call(rbind, results_modal_cont_centrality_mlvar)

results_modal_cont_centrality_mlvar <- as.data.frame(results_modal_cont_centrality_mlvar)

results_modal_cont_centrality_mlvar$Participant <- rownames(results_modal_cont_centrality_mlvar)

rownames(results_modal_cont_centrality_mlvar) <- NULL

results_modal_cont_centrality_mlvar <- results_modal_cont_centrality_mlvar %>% relocate(Participant, .before = V1)

colnames(results_modal_cont_centrality_mlvar) <- c("Participant", "bad_mcc", "control_mcc", "energy_mcc", "focus_mcc", 
                                                   "fun_mcc", "interest_mcc", "movement_mcc", "sad_mcc")

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
# Computing expected influence for all participants for the ML-VAR 7-node network with control ----
# ---------------------------------------------------------------------------- #

results_exp_inf_mlvar_control <- lapply(adjacency_matrices_mlvar_control,
                                        function(matrix) expectedInf(matrix, step = c("both", 1, 2), directed = TRUE))

#transform this list into a data frame

expected_influence_mlvar_control <- data.frame()

for (participant in names(results_exp_inf_mlvar_control)) {
  participant_data <- results_exp_inf_mlvar_control[[participant]]

  for (step in names(participant_data)) {
    step_data <- participant_data[[step]]

    temp_df <- data.frame(
      Participant = participant,
      Step = step,
      t(as.data.frame(step_data))
    )

    expected_influence_mlvar_control <- rbind(expected_influence_mlvar_control, temp_df)
  }
}

rownames(expected_influence_mlvar_control) <- NULL

expected_influence_mlvar_control_wide <- expected_influence_mlvar_control %>%
  pivot_wider(
    names_from = Step, 
    values_from = c(bad, control, energy, focus, interest, movement, sad),
    names_sep = "_"
  ) 

expected_influence_mlvar_control_wide <- expected_influence_mlvar_control_wide [, c("Participant", 
                                                                                    "bad_step1", "control_step1", "energy_step1", "focus_step1", "interest_step1", "movement_step1", "sad_step1", 
                                                                                    "bad_step2", "control_step2", "energy_step2", "focus_step2", "interest_step2", "movement_step2", "sad_step2")]

expected_influence_just_mlvar_control <- expected_influence_mlvar_control_wide [, c("Participant", "control_step1", "control_step2")]

colnames(expected_influence_just_mlvar_control) <- c("Participant", "control_step1_7", "control_step2_7")

# ---------------------------------------------------------------------------- #
# Computing average controllability centrality for the 7-node ML-VAR network with control ----
# ---------------------------------------------------------------------------- #

#ave_mlvar_control_centrality function

#Computing average control centrality for all participants and transforming the resulting list into a dataframe

results_ave_cont_centrality_mlvar_control <- lapply(adjacency_matrices_mlvar_control, function(matrix) ave_control_centrality(matrix))

results_ave_cont_centrality_mlvar_control <- do.call(rbind, results_ave_cont_centrality_mlvar_control)

results_ave_cont_centrality_mlvar_control <- as.data.frame(results_ave_cont_centrality_mlvar_control)

results_ave_cont_centrality_mlvar_control$Participant <- rownames(results_ave_cont_centrality_mlvar_control)

rownames(results_ave_cont_centrality_mlvar_control) <- NULL

results_ave_cont_centrality_mlvar_control <- results_ave_cont_centrality_mlvar_control %>% relocate(Participant, .before = V1)

colnames(results_ave_cont_centrality_mlvar_control) <- c("Participant", "bad_acc", "control_acc", "energy_acc", "focus_acc", 
                                                         "interest_acc", "movement_acc", "sad_acc")

results_ave_cont_centrality_just_mlvar_control <- results_ave_cont_centrality_mlvar_control [, c("Participant", "control_acc")]

colnames(results_ave_cont_centrality_just_mlvar_control) <- c("Participant", "control_acc_7")

# ---------------------------------------------------------------------------- #
# Computing modal controllability centrality for the ML-VAR 7-node network with control ----
# ---------------------------------------------------------------------------- #

results_modal_cont_centrality_mlvar_control <- lapply(adjacency_matrices_mlvar_control, function(matrix) modal_control_centrality(matrix))

results_modal_cont_centrality_mlvar_control <- do.call(rbind, results_modal_cont_centrality_mlvar_control)

results_modal_cont_centrality_mlvar_control <- as.data.frame(results_modal_cont_centrality_mlvar_control)

results_modal_cont_centrality_mlvar_control$Participant <- rownames(results_modal_cont_centrality_mlvar_control)

rownames(results_modal_cont_centrality_mlvar_control) <- NULL

results_modal_cont_centrality_mlvar_control <- results_modal_cont_centrality_mlvar_control %>% relocate(Participant, .before = V1)

colnames(results_modal_cont_centrality_mlvar_control) <- c("Participant", "bad_mcc", "control_mcc", "energy_mcc", "focus_mcc", 
                                                           "interest_mcc", "movement_mcc", "sad_mcc")

results_modal_cont_centrality_just_mlvar_control <- results_modal_cont_centrality_mlvar_control[, c("Participant", "control_mcc")]

colnames(results_modal_cont_centrality_just_mlvar_control) <- c("Participant", "control_mcc_7")

# ---------------------------------------------------------------------------- #
# Computing expected influence for all participants for the ML-VAR 7-node network with fun ----
# ---------------------------------------------------------------------------- #

results_exp_inf_mlvar_fun <- lapply(adjacency_matrices_mlvar_fun,
                                        function(matrix) expectedInf(matrix, step = c("both", 1, 2), directed = TRUE))

#transform this list into a data frame

expected_influence_mlvar_fun <- data.frame()

for (participant in names(results_exp_inf_mlvar_fun)) {
  participant_data <- results_exp_inf_mlvar_fun[[participant]]
  
  for (step in names(participant_data)) {
    step_data <- participant_data[[step]]
    
    temp_df <- data.frame(
      Participant = participant,
      Step = step,
      t(as.data.frame(step_data))
    )
    
    expected_influence_mlvar_fun <- rbind(expected_influence_mlvar_fun, temp_df)
  }
}

rownames(expected_influence_mlvar_fun) <- NULL

expected_influence_mlvar_fun_wide <- expected_influence_mlvar_fun %>%
  pivot_wider(
    names_from = Step, 
    values_from = c(bad, energy, focus, fun, interest, movement, sad),
    names_sep = "_"
  ) 

expected_influence_mlvar_fun_wide <- expected_influence_mlvar_fun_wide [, c("Participant", 
                                                                            "bad_step1", "energy_step1", "focus_step1", "fun_step1", "interest_step1", "movement_step1", "sad_step1", 
                                                                            "bad_step2", "energy_step2", "focus_step2", "fun_step2", "interest_step2", "movement_step2", "sad_step2")]

expected_influence_just_mlvar_fun <- expected_influence_mlvar_fun_wide [, c("Participant", "fun_step1", "fun_step2")]

colnames(expected_influence_just_mlvar_fun) <- c("Participant", "fun_step1_7", "fun_step2_7")

# ---------------------------------------------------------------------------- #
# Computing average controllability centrality for the ML-VAR 7-node network with fun ----
# ---------------------------------------------------------------------------- #

#ave_mlvar_fun_centrality function

#Computing average control centrality for all participants and transforming the resulting list into a dataframe

results_ave_cont_centrality_mlvar_fun <- lapply(adjacency_matrices_mlvar_fun, function(matrix) ave_control_centrality(matrix))

results_ave_cont_centrality_mlvar_fun <- do.call(rbind, results_ave_cont_centrality_mlvar_fun)

results_ave_cont_centrality_mlvar_fun <- as.data.frame(results_ave_cont_centrality_mlvar_fun)

results_ave_cont_centrality_mlvar_fun$Participant <- rownames(results_ave_cont_centrality_mlvar_fun)

rownames(results_ave_cont_centrality_mlvar_fun) <- NULL

results_ave_cont_centrality_mlvar_fun <- results_ave_cont_centrality_mlvar_fun %>% relocate(Participant, .before = V1)

colnames(results_ave_cont_centrality_mlvar_fun) <- c("Participant", "bad_acc", "energy_acc", "focus_acc", "fun_acc",
                                                     "interest_acc", "movement_acc", "sad_acc")

results_ave_cont_centrality_just_mlvar_fun <- results_ave_cont_centrality_mlvar_fun [, c("Participant", "fun_acc")]

colnames(results_ave_cont_centrality_just_mlvar_fun) <- c("Participant", "fun_acc_7")

# ---------------------------------------------------------------------------- #
# Computing modal controllability centrality for the ML-VAR 7-node network with fun ----
# ---------------------------------------------------------------------------- #

results_modal_cont_centrality_mlvar_fun <- lapply(adjacency_matrices_mlvar_fun, function(matrix) modal_control_centrality(matrix))

results_modal_cont_centrality_mlvar_fun <- do.call(rbind, results_modal_cont_centrality_mlvar_fun)

results_modal_cont_centrality_mlvar_fun <- as.data.frame(results_modal_cont_centrality_mlvar_fun)

results_modal_cont_centrality_mlvar_fun$Participant <- rownames(results_modal_cont_centrality_mlvar_fun)

rownames(results_modal_cont_centrality_mlvar_fun) <- NULL

results_modal_cont_centrality_mlvar_fun <- results_modal_cont_centrality_mlvar_fun %>% relocate(Participant, .before = V1)

colnames(results_modal_cont_centrality_mlvar_fun) <- c("Participant", "bad_mcc", "energy_mcc", "focus_mcc", "fun_mcc",
                                                       "interest_mcc", "movement_mcc", "sad_mcc")

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