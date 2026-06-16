library(tidyverse)
library(tidybayes)
library(rstan)
seed <- 404806
set.seed(seed)
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
n_draws <- 5
file_path <- "/C_auris_OD_output/fake_data_check"
#
# Get job idx
#
full_job_idx_list <- expand.grid(Draw = 1:n_draws, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up prev job list
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
full_job_idx_list %>%
  filter(Model == job_idx) %>%
  select(-Model) -> job_idx_list
#
# get model
#
model <- model_list[job_idx]
#
# read in model fits
#
fit_idx <- which(full_job_idx_list$Model==job_idx)
stan_list <- lapply(fit_idx, function(x) readRDS(file = paste(storage_loc, file_path, "/", model, "_", x, "_fit.Rda", sep="")))
#
# combine stan fits from chains but using the same fake_data into a single stan fit
#
stan_fit <- lapply(1:n_draws, function(x) sflist2stanfit(stan_list[job_idx_list$Draw==x])) 
#
# save merged posterior tibbles
#
lapply(1:n_draws, function(i) saveRDS(gather_draws(stan_fit[[i]], cbind(gamma, kappa, epsilon)[k]) , paste("output/fake_data_check/", model, "_posterior_params_", i, ".Rda", sep = "")))  

