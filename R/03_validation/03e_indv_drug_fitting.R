library(tidyverse)
library(rstan)
library(tidybayes)
seed <- 404806
set.seed(seed)
model_list <- c("edwards_indv_drug", 
                "gompertz_indv_drug",
                "edwards_indv_drug_epsilon0",
                "edwards_indv_drug_gamma0", 
                "edwards_indv_drug_kappa0",
                "edwards_indv_drug_decay", 
                "edwards_indv_drug_hs_reg")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
no_of_folds <- 3
file_path <- "/C_auris_OD_output/indv_drug"
iter <- 2000 
#
# Create idx list
#
job_idx_list <- expand.grid(Fold = 1:no_of_folds, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
#
# Get job index
#
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
#
# Specify labelled data
#
fold_df <- readRDS(file=paste(storage_loc, file_path, "/data/fold_data_", job_idx, ".Rda", sep=""))
#
# Specify model
#
model_idx <- job_idx_list$Model[job_idx]
model <- model_list[model_idx]
stanfile <- paste('models/', model, '.stan', sep = "")
stanmodel <- stan_model(stanfile)
#
# create stan data
#
fold_df %>%
  select(drug_conds, drug1_name, drug2_name) %>%
  unique() %>%
  arrange(drug_conds) %>% 
  model.matrix(~drug1_name+drug2_name, data=.) -> X  # matrix of covariates for added drugs
colnames(X) <- NULL
# create matrix of indicators of drug conditions x drugs (AFG, MGX and 5FC)
single_drug_matrix <- cbind(rep(1, dim(X)[1]), X[,2:length(unique(c(fold_df$drug1_name, fold_df$drug2_name)))] + X[,(length(unique(c(fold_df$drug1_name, fold_df$drug2_name)))+1):(length(unique(c(fold_df$drug1_name, fold_df$drug2_name)))*2-1)])

stan_data <- list(D = length(unique(c(fold_df$drug1_name, fold_df$drug2_name))),
                  D_c = length(unique(fold_df$drug_conds)),  
                  X = single_drug_matrix,
                    
                  N_train = length(which(fold_df$training)),
                  T_train = length(unique(filter(fold_df, training)$time)),
                  train_time = sort(unique(filter(fold_df, training)$time)),
                  time_train_idx = match(filter(fold_df, training)$time, sort(unique(filter(fold_df, training)$time))), 
                  max_train_time = max(sort(unique(filter(fold_df, training)$time))),
                  N_b_train = length(which(filter(fold_df, training)$blanks)),
                  blank_train_idx = which(filter(fold_df, training)$blanks),
                  N_auris_train = length(which(!filter(fold_df, training)$blanks)),
                  auris_train_idx = which(!filter(fold_df, training)$blanks),
                  y_obs_train = filter(fold_df, training)$OD,
                  drug_cond_train_idx = as.numeric(filter(fold_df, training)$drug_conds),
                  
                  N_test = length(which(!fold_df$training)),
                  T_test = length(unique(filter(fold_df, testing)$time)),
                  test_time = as.array(sort(unique(filter(fold_df, testing)$time))),
                  time_test_idx = match(filter(fold_df, testing)$time, sort(unique(filter(fold_df, testing)$time))),
                  N_b_test = length(which(filter(fold_df, testing)$blanks)),
                  blank_test_idx = which(filter(fold_df, testing)$blanks),
                  N_auris_test = length(which(!filter(fold_df, testing)$blanks)),
                  auris_test_idx = which(!filter(fold_df, testing)$blanks),
                  y_obs_test = filter(fold_df, !training)$OD,
                  drug_cond_test_idx = as.numeric(filter(fold_df, testing)$drug_conds),
                  
                  y_train_idx = which(fold_df$training),
                  y_test_idx = which(!fold_df$training),
                  include_likelihood = 1)

# Fit Model --------------------------------------------------
chosen_chain <- job_idx_list$Chains[job_idx]
stan_fit <- sampling(stanmodel, data = stan_data, chains = 1, iter = iter, seed = seed, chain_id = chosen_chain, cores = 1, control = list(max_treedepth = 12))
if (!dir.exists(file.path(storage_loc, file_path))) {
  dir.create(file.path(storage_loc, file_path))
}
saveRDS(stan_fit, file=paste(storage_loc, file_path, "/", model, "_", job_idx, "_fit.Rda", sep=""))  #stores model
