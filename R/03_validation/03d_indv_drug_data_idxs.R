library(tidyverse)
seed <- 404806
set.seed(seed)
storage_loc <- Sys.getenv("EPHEMERAL")  #store fold data in TMP_DIR
model_list <- c("edwards_indv_drug", 
                "gompertz_indv_drug",
                "edwards_indv_drug_epsilon0",
                "edwards_indv_drug_gamma0", 
                "edwards_indv_drug_kappa0",
                "edwards_indv_drug_decay", 
                "edwards_indv_drug_hs_reg")
#
# Specify script options
#
no_of_chains <- 4
# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD.Rda") %>%
  drop_na(drug_conds) %>%  # drop wells with nothing in them
  mutate(drug_conds = factor(drug_conds, levels = c("2X RPMI + dH2O", 
                                                    "RPMI", 
                                                    "8 AFG", 
                                                    "0.03 MGX", 
                                                    "0.25 5FC", 
                                                    "8 AFG +\r\n0.03 MGX", 
                                                    "8 AFG +\r\n0.25 5FC"))) %>%
  filter(!drug_conds%in%c("8 AFG +\r\n0.03 MGX", "8 AFG +\r\n0.25 5FC")) %>%  # remove combination drug conditions
  select(-strain) %>%
  separate(drug_name, c("drug1_name", "drug2_name"), "_", remove=TRUE) %>%
  pivot_longer(c(drug1_name, drug2_name), names_to="drug", values_to = "drug_name") %>% 
  mutate(drug_name = factor(drug_name, levels = c(NA, "AFG", "MGX", "5FC"), labels = c("None", "AFG", "MGX", "5FC"), exclude=NULL)) %>%
  pivot_wider(names_from = "drug", values_from = "drug_name") -> data
#
# Create idx list
# one well from each drug condition is used as held-out testing data for each fold
#
data %>% 
  filter(!blanks) %>% 
  group_by(drug_conds) %>% 
  select(well) %>% 
  unique() %>% 
  mutate(fold = as.numeric(factor(well))) %>%  # each rep in drug conditions given fold
  mutate(fold = sample(fold)) %>%  # shuffle allocated folds within each group
  ungroup() %>%
  select(well, fold) -> test_wells 
no_of_folds <- max(test_wells$fold)
job_idx_list <- expand.grid(Fold = 1:no_of_folds, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
labelled_data <- lapply(1:nrow(job_idx_list), function(x){data %>%
                                                            mutate(testing = well%in%filter(test_wells, fold==job_idx_list$Fold[x])$well) %>% 
                                                            mutate(training = !testing)})
file_path <- paste(storage_loc, "/C_auris_OD_output/indv_drug/data", sep="")
if (!dir.exists(file.path(file_path))) {
  dir.create(file.path(file_path))
}
lapply(1:nrow(job_idx_list), function(x){saveRDS(labelled_data[[x]], file=paste(file_path, "/fold_data_", x, ".Rda", sep=""))})  #stores fold data
