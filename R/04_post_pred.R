#
# Script that runs a posterior predictive check (to be visually assessed) 
#
rm(list = ls())
library(rstan)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
seed <- 404806
set.seed(seed)

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
  pivot_wider(names_from = "drug", values_from = "drug_name") %>%
  add_column(training = TRUE) %>%
  mutate(testing = !training) -> data

# Specify and fit model ---------------------------------------------------
model <- "edwards_drug_hs_reg_full_data"   # switch with "gompertz_drug_hs_reg_full_data" to get the gompertz fit
stanfile <- paste("models/", model, ".stan", sep = "")
no_of_chains <- 4
iter <- 2000
#
# sets up stan data
#
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
#
#  sets up stan data in same order as in edwards_drug_hs_reg_full_data and gompertz_drug_hs_reg_full_data models
#
stan_data <- list(D = length(unique(c(data$drug1_name, data$drug2_name))),
                  D_i = dim(interaction_drug_matrix)[2],
                  D_c = length(unique(data$drug_conds)),  
                  X = cbind(single_drug_matrix, interaction_drug_matrix),  # combine binary matrices for drug and drug interactions
                  
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
                  
                  include_likelihood = 1)
# Fits model
fit <- stan(file = stanfile,
            data = stan_data,
            seed = seed,
            chains = no_of_chains,
            iter = iter,
            control=list(max_treedepth=12))

#
# save model
#
saveRDS(fit, file=paste("output/", model, "_full_fit.Rda", sep=""))  #stores model
