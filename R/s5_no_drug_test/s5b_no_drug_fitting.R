library(tidyverse)
library(rstan)
library(tidybayes)
seed <- 404806
set.seed(seed)
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
no_of_trials <- 3
file_path <- "/C_auris_OD_output/no_drug_test"
iter <- 2000 
#
# Create idx list
#
job_idx_list <- expand.grid(Trial = 1:no_of_trials, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
#
# Get job index
#
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
#
# Specify labelled data
#
fake_df <- readRDS(file=paste(storage_loc, file_path, "/data/fake_data_", job_idx, ".Rda", sep=""))
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
fake_df %>%
  select(drug_conds, drug_name) %>%
  unique() %>%
  arrange(drug_conds) %>% 
  model.matrix(~drug_name, data=.) -> X  # matrix of covariates for drug 
colnames(X) <- NULL

stan_data <- list(D = length(unique(fake_df$drug_name)),
                  D_i = 0, # only considering one test drug here
                  D_c = length(unique(fake_df$drug_conds)),  
                  X = X,   

                  N = nrow(fake_df),
                  T = length(unique(fake_df$time)),
                  time = sort(unique(fake_df$time)),
                  time_idx = match(fake_df$time, sort(unique(fake_df$time))),  
                  max_time = max(sort(unique(fake_df$time))),
                  N_b = length(which(fake_df$blanks)),
                  blank_idx = which(fake_df$blanks),
                  N_auris = length(which(!fake_df$blanks)),
                  auris_idx = which(!fake_df$blanks),
                  y_obs = fake_df$OD,
                  drug_cond_idx = as.numeric(fake_df$drug_conds),

                  include_likelihood = 1)

# Fit Model --------------------------------------------------
chosen_chain <- job_idx_list$Chains[job_idx]
stan_fit <- sampling(stanmodel, data = stan_data, chains = 1, iter = iter, seed = seed, chain_id=chosen_chain, cores = 1, control = list(max_treedepth = 12))
if (!dir.exists(file.path(storage_loc, file_path))) {
  dir.create(file.path(storage_loc, file_path))
}
saveRDS(stan_fit, file=paste(storage_loc, file_path, "/", model, "_", job_idx, "_fit.Rda", sep=""))  #stores model

