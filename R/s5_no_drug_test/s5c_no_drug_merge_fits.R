library(tidyverse)
library(rstan)
library(tidybayes)
seed <- 404806
set.seed(seed)
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
#
# Read in data set
#
readRDS(file = "data/OD.Rda") %>%
  filter(drug_conds%in%c("RPMI", "2X RPMI + dH2O")) %>%  #filter for drug free data 
  select(-c(drug_conc, strain)) -> data
#
# Get fake data
#
non_blank_wells <- unique(filter(data, !blanks)$well)
no_of_trials <- length(non_blank_wells)
labelled_data <- lapply(1:no_of_trials, function(x){data %>%
                                                      mutate(trial = x) %>% 
                                                      mutate(drug_conds = case_when(well == non_blank_wells[x] ~ "test_drug", .default = drug_conds)) %>%
                                                      mutate(drug_name = case_when(drug_conds == "test_drug" ~ "test", .default = "None")) %>%
                                                      mutate(drug_conds = factor(drug_conds, levels = c("2X RPMI + dH2O", "RPMI", "test_drug")))})
file_path <- "/C_auris_OD_output/no_drug_test"
#
# Get job idx
#
full_job_idx_list <- expand.grid(Trial = 1:no_of_trials, Chains = 1:no_of_chains, Model = 1:length(model_list)) #prev job list used in s5a_
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
full_job_idx_list %>%
  filter(Model == job_idx) %>%
  select (-Model) -> job_idx_list
#
# get model
#
model <- model_list[job_idx]
#
# read in model training fits
#
fit_idx <- which(full_job_idx_list$Model==job_idx)
stan_list <- lapply(fit_idx, function(x) readRDS(file = paste(storage_loc, file_path, "/", model, "_", x, "_fit.Rda", sep="")))
#
# merge different chains from same fake data fit
#
stan_fit <- lapply(1:no_of_trials, function(x) sflist2stanfit(stan_list[job_idx_list$Trial==x])) 
#
# get and save drug posterior params
#
lapply(1:no_of_trials, function(i) saveRDS(gather_draws(stan_fit[[i]], cbind(gamma, kappa, epsilon)[k]) , paste("output/no_drug_test/", model, "_posterior_params_", i, ".Rda", sep = "")))  
