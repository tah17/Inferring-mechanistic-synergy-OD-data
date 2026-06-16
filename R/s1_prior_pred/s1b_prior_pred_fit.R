library(rstan)
library(tidyverse)
library(tidybayes)
seed <- 404806
set.seed(seed)
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
storage_loc <- Sys.getenv("EPHEMERAL")  #store large stan fits in TMP_DIR
no_of_chains <- 4
n_draws <- 5
file_path <- "/C_auris_OD_output/fake_data_check/"
iter <- 2000
#
# Create idx list
#
job_idx_list <- expand.grid(Draw = 1:n_draws, Chains = 1:no_of_chains, Model = 1:length(model_list))  #sets up job list
#
# Get job index
#
job_idx <- as.integer(Sys.getenv("PBS_ARRAY_INDEX"))
#
# Read in fake data
#
fake_data <- readRDS(file=paste(storage_loc, file_path, "/data/fake_data_", job_idx, ".Rda", sep=""))
#
# Specify model
#
model_idx <- job_idx_list$Model[job_idx]
model <- model_list[model_idx]
stanfile <- paste('models/', model, '.stan', sep = "")
stanmodel <- stan_model(stanfile)

# Specify fake data ---------------------------------------------------
fake_data %>%
  select(drug_conds, drug1_name, drug2_name) %>%
  unique() %>%
  arrange(drug_conds) %>% 
  model.matrix(~drug1_name+drug2_name, data=.) -> X  # matrix of covariates for drug 1 and drug 2
colnames(X) <- NULL
# create matrix of indicators of drug conditions x drugs (AFG, MGX and 5FC)
single_drug_matrix <- cbind(rep(1, dim(X)[1]), X[,2:length(unique(c(fake_data$drug1_name, fake_data$drug2_name)))] + X[,(length(unique(c(fake_data$drug1_name, fake_data$drug2_name)))+1):(length(unique(c(fake_data$drug1_name, fake_data$drug2_name)))*2-1)])
# create matrix of indicators of drug conditions x drug interactions (AFG:MGX and AFG:5FC)
interaction_drug_matrix <- sapply(combn(2:length(unique(c(fake_data$drug1_name, fake_data$drug2_name))), 2, simplify=FALSE), function(x) single_drug_matrix[, x[1]]*single_drug_matrix[, x[2]])
interaction_drug_matrix <- interaction_drug_matrix[, -(length(unique(c(fake_data$drug1_name, fake_data$drug2_name)))-1)]

stan_data <- list(D = length(unique(c(fake_data$drug1_name, fake_data$drug2_name))),
                  D_i = dim(interaction_drug_matrix)[2],
                  D_c = length(unique(fake_data$drug_conds)),  
                  X = cbind(single_drug_matrix, interaction_drug_matrix),   # combine binary matrices for drug and drug interactions

                  N = nrow(fake_data),
                  T = length(unique(fake_data$time)),
                  time = sort(unique(fake_data$time)),
                  time_idx = match(fake_data$time, sort(unique(fake_data$time))),  
                  max_time = max(sort(unique(fake_data$time))),
                  N_b = length(which(fake_data$blanks)),
                  blank_idx = which(fake_data$blanks),
                  N_auris = length(which(!fake_data$blanks)),
                  auris_idx = which(!fake_data$blanks),
                  y_obs = fake_data$OD,
                  drug_cond_idx = as.numeric(fake_data$drug_conds),
                  
                  include_likelihood = 1)

# Fit Model --------------------------------------------------
chosen_chain <- job_idx_list$Chains[job_idx]
# fit model to its own fake data
stan_fit <- sampling(stanmodel, data = stan_data, chains = 1, iter = iter, seed = seed, chain_id=chosen_chain, cores = 1, control = list(max_treedepth = 12, adapt_delta=0.99))
if (!dir.exists(file.path(storage_loc, file_path))) {
  dir.create(file.path(storage_loc, file_path))
}
saveRDS(stan_fit, file=paste(storage_loc, file_path, "/", model, "_", job_idx, "_fit.Rda", sep=""))  #stores model

