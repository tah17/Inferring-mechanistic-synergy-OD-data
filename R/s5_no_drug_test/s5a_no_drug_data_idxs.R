library(tidyverse)
seed <- 404806
set.seed(seed)
storage_loc <- Sys.getenv("EPHEMERAL")  #store fold data in TMP_DIR
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
#
# Specify script options
#
no_of_chains <- 4
# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD.Rda") %>%
  filter(drug_conds%in%c("RPMI", "2X RPMI + dH2O")) %>%  #filter for drug free data 
  select(-c(drug_conc, strain)) -> data
#
# Create idx list
#
non_blank_wells <- unique(filter(data, !blanks)$well)
job_idx_list <- expand.grid(Trial = 1:length(non_blank_wells), Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
labelled_data <- lapply(1:nrow(job_idx_list), function(x){data %>%
                                                              mutate(drug_conds = case_when(well == non_blank_wells[job_idx_list$Trial[x]] ~ "test_drug", .default = drug_conds)) %>%
                                                              mutate(drug_name = case_when(drug_conds == "test_drug" ~ "test", .default = "None")) %>%
                                                              mutate(drug_conds = factor(drug_conds, levels = c("2X RPMI + dH2O", "RPMI", "test_drug")))})
file_path <- paste(storage_loc, "/C_auris_OD_output/no_drug_test/data", sep="")

if (!dir.exists(file.path(file_path))) {
  dir.create(file.path(file_path))
}
lapply(1:nrow(job_idx_list), function(x){saveRDS(labelled_data[[x]], file=paste(file_path, "/fake_data_", x, ".Rda", sep=""))})  #stores fake data

