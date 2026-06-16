library(tidyverse)
seed <- 404806
set.seed(seed)
storage_loc <- Sys.getenv("EPHEMERAL")  #store fold data in TMP_DIR
model_list <- c("edwards_OD_direct", 
                "edwards", 
                "exponential", 
                "logistic", 
                "gompertz", 
                "gompertz_OD_direct")
#
# Specify script options
#
no_of_chains <- 4
# Read in Data ------------------------------------------------------------
readRDS(file = "data/OD.Rda") %>%
  filter(drug_conds%in%c("RPMI", "2X RPMI + dH2O")) %>%  #filter for drug free data 
  select(-c(drug_name, drug_conc, strain, drug_conds)) -> data
#
# Create idx list
#
no_of_folds <- length(unique(filter(data, !blanks)$well))
test_wells <- unique(filter(data, !blanks)$well)
job_idx_list <- expand.grid(Fold = 1:no_of_folds, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
labelled_data <- lapply(1:nrow(job_idx_list), function(x){data %>%
                                                              mutate(testing = well == test_wells[job_idx_list$Fold[x]]) %>% 
                                                              mutate(training = !testing)})
file_path <- paste(storage_loc, "/C_auris_OD_output/drug_free/data", sep="")

if (!dir.exists(file.path(file_path))) {
  dir.create(file.path(file_path))
}
lapply(1:nrow(job_idx_list), function(x){saveRDS(labelled_data[[x]], file=paste(file_path, "/fold_data_", x, ".Rda", sep=""))})  #stores fold data
