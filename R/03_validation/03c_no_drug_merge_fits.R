library(tidyverse)
library(rstan)
seed <- 404806
set.seed(seed)
source("R/03_validation/functions.R")
model_list <- c("edwards_OD_direct", 
                "edwards", 
                "exponential", 
                "logistic", 
                "gompertz", 
                "gompertz_OD_direct")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
#
# Read in data set
#
readRDS(file = "data/OD.Rda") %>%
  filter(drug_conds%in%c("RPMI", "2X RPMI + dH2O")) %>%
  select(-c(drug_name, drug_conc, strain, drug_conds)) -> data
#
# Get labelled data
#
no_of_folds <- length(unique(filter(data, !blanks)$well))
test_wells <- unique(filter(data, !blanks)$well)
labelled_data <- lapply(1:no_of_folds, function(x){data %>% 
                                                    mutate(fold = x) %>% 
                                                    mutate(testing = well == test_wells[x]) %>% 
                                                    mutate(training = !testing)})
file_path <- "/C_auris_OD_output/drug_free"
#
# Get job idx
#
full_job_idx_list <- expand.grid(Fold = 1:no_of_folds, Chains = 1:no_of_chains, Model = 1:length(model_list))  #prev job list used in 03a_
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
full_job_idx_list %>%
  filter(Model == job_idx) %>%
  select(-Model) -> job_idx_list
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
# merge different chains from same fold
#
# combine stan fits from chains but using the same cv fold into a single stan fit
stan_fit <- lapply(1:no_of_folds, function(x) sflist2stanfit(stan_list[job_idx_list$Fold==x])) 
#
# save merged stan fit objects
#
lapply(1:no_of_folds, function(i) saveRDS(stan_fit[[i]], paste(storage_loc, file_path, "/fold_", i, "_", model, "_train.rds", sep = "")))  
#
# calculate training + testing errors and save
#
fit_list <- lapply(1:no_of_folds, function(x) get_metrics(labelled_data[[x]], stan_fit[[x]]) %>% mutate(fold = x)) 
saveRDS(do.call(rbind, fit_list), paste("output/drug_free/", model, "_fit_stats.rds", sep = ""))

