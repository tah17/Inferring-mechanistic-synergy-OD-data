library(tidyverse)
seed <- 404806
set.seed(seed)
storage_loc <- Sys.getenv("EPHEMERAL")  #store fold data in TMP_DIR
#
# Specify script options
#
source("R/s1_prior_pred/functions.R")
iter <- 2000
n_draws <- 5   # number of times to do fake data check per model
no_of_chains <- 4
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
job_idx_list <- expand.grid(Draw = 1:n_draws, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
# Read in data ------------------------------------------------------------
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

# Specify and fit model ---------------------------------------------------
data %>%
  select(drug_conds, drug1_name, drug2_name) %>%
  unique() %>%
  arrange(drug_conds) %>% 
  model.matrix(~drug1_name+drug2_name, data=.) -> X  # matrix of covariates for drug 1 and drug 2
colnames(X) <- NULL
# create matrix of indicators of drug conditions x drugs (AFG, MGX and 5FC)
single_drug_matrix <- cbind(rep(1, dim(X)[1]), X[,2:length(unique(c(data$drug1_name, data$drug2_name)))] + X[,(length(unique(c(data$drug1_name, data$drug2_name)))+1):(length(unique(c(data$drug1_name, data$drug2_name)))*2-1)])
# create matrix of indicators of drug conditions x drug interactions (AFG:MGX and AFG:5FC)
interaction_drug_matrix <- sapply(combn(2:length(unique(c(data$drug1_name, data$drug2_name))), 2, simplify=FALSE), function(x) single_drug_matrix[, x[1]]*single_drug_matrix[, x[2]])
interaction_drug_matrix <- interaction_drug_matrix[, -(length(unique(c(data$drug1_name, data$drug2_name)))-1)]

stan_data <- list(D = length(unique(c(data$drug1_name, data$drug2_name))),
                  D_i = dim(interaction_drug_matrix)[2],
                  D_c = length(unique(data$drug_conds)),  
                  X = cbind(single_drug_matrix, interaction_drug_matrix),   # combine binary matrices for drug and drug interactions

                  N = nrow(data),
                  T = length(unique(data$time)),
                  time = sort(unique(data$time)),
                  time_idx = match(data$time, sort(unique(data$time))),  
                  max_time = max(sort(unique(data$time))),
                  N_b = length(which(data$blanks)),
                  blank_idx = which(data$blanks),
                  N_auris = length(which(!data$blanks)),
                  auris_idx = which(!data$blanks),
                  y_obs = data$OD,
                  drug_cond_idx = as.numeric(data$drug_conds),
                  
                  include_likelihood = 0)  # include_likelihood set to 0 to sample from prior
# Get prior, fake data and drug-action parameters used to generate fake data
prior_samples <- vector("list", length(model_list)) 
fake_data <- vector("list", length(model_list))
fake_data_params <- vector("list", length(model_list))
for (i in seq_along(model_list)) {
  stanfile <- paste("models/", model_list[i], ".stan", sep = "")
  stanmodel <- stan_model(stanfile)
  prior_samples[[i]] <- get_prior_samples(stanmodel, stan_data, no_of_chains, iter, seed, data) 
  fake_data_draws <- lapply(1:n_draws, function(x){get_fake_data(prior_samples[[i]], data)})
  fake_data[[i]] <- lapply(fake_data_draws, function(x) x[[1]])
  fake_data_params[[i]] <- lapply(fake_data_draws, function(x) x[[2]])
}
file_path <- paste(storage_loc, "/C_auris_OD_output/fake_data_check/data", sep="")
if (!dir.exists(file.path(file_path))) {
  dir.create(file.path(file_path))
}
filter(job_idx_list, Chains==1) %>%
  select(-Chains) -> draw_model_list
lapply(1:length(model_list), function(x){saveRDS(prior_samples[[x]]$drug_priors, file=paste("output/fake_data_check/drug_priors_", model_list[x], ".Rda", sep=""))})  #stores priors
lapply(1:nrow(job_idx_list), function(x){saveRDS(fake_data[[job_idx_list$Model[x]]][[job_idx_list$Draw[x]]], file=paste(file_path, "/fake_data_", x, ".Rda", sep=""))})  #stores fake data
lapply(1:nrow(draw_model_list), function(x){saveRDS(fake_data_params[[draw_model_list$Model[x]]][[draw_model_list$Draw[x]]], file=paste("output/fake_data_check/fake_data_params_", x, ".Rda", sep=""))})  #stores fake data params

