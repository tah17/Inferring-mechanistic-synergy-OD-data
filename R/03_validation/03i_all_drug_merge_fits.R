library(tidyverse)
library(rstan)
seed <- 404806
set.seed(seed)
source("R/03_validation/functions.R")
model_list <- c("edwards_drug_no_syn_hs_reg",
                "edwards_drug_hs_reg", 
                "gompertz_drug_hs_reg",
                "gompertz_drug_no_syn_hs_reg")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
#
# Read in data set
#
readRDS(file = "data/OD.Rda") %>%
  drop_na(drug_conds) %>%  # drop wells with nothing in them
  mutate(drug_conds = factor(drug_conds, levels = c("2X RPMI + dH2O", 
                                                    "RPMI", 
                                                    "8 AFG", 
                                                    "0.03 MGX", 
                                                    "0.25 5FC", 
                                                    "8 AFG +\r\n0.03 MGX", 
                                                    "8 AFG +\r\n0.25 5FC"))) %>%
  select(-strain) %>%
  separate(drug_name, c("drug1_name", "drug2_name"), "_", remove=TRUE) %>%
  pivot_longer(c(drug1_name, drug2_name), names_to="drug", values_to = "drug_name") %>% 
  mutate(drug_name = factor(drug_name, levels = c(NA, "AFG", "MGX", "5FC"), labels = c("None", "AFG", "MGX", "5FC"), exclude=NULL)) %>%
  pivot_wider(names_from = "drug", values_from = "drug_name") -> data
#
# Get labelled data
#
data %>% 
  filter(!blanks) %>% 
  group_by(drug_conds) %>% 
  select(well) %>% 
  unique() %>% 
  mutate(fold = as.numeric(factor(well))) %>%
  mutate(fold = sample(fold)) %>%
  ungroup() %>%
  select(well, fold) -> test_wells 
no_of_folds <- max(test_wells$fold)
labelled_data <- lapply(1:no_of_folds, function(x){data %>% 
                                                    mutate(fold = x) %>% 
                                                    mutate(testing = well%in%filter(test_wells, fold==x)$well) %>% 
                                                    mutate(training = !testing)})
file_path <- "/C_auris_OD_output/all_drug"
#
# Get job idx
#
full_job_idx_list <- expand.grid(Fold = 1:no_of_folds, Chains = 1:no_of_chains, Model = 1:length(model_list))  #prev job list used in 03h_
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
saveRDS(do.call(rbind, fit_list), paste("output/all_drug/", model, "_fit_stats.rds", sep = ""))
