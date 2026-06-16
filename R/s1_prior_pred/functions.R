library(rstan)
library(tidyverse)
library(tidybayes)

#' get_fake_data
#' 
#' This function takes the output of the get_prior_samples function and returns some fake data and the drug parameters that generated the fake data 
#' 
#' @param prior_samples List. List returned get_prior_samples() below.
#' @param data Tibble. Tibble of real data.
#' 
#' @return List of Tibble of fake data and Tibble of true parameters used to generate the fake data.
#' 
get_fake_data <- function(prior_samples, data) {
  # Draws one sample from the prior predictive distribution that is not flat
  draw <- sample(prior_samples$viable_draws, 1)
  gather_draws(prior_samples$prior, y_tot[1]) %>%
    filter(.draw == draw) -> y_rep_draws
  y_rep_draw <- y_rep_draws$.value[[1]]
  
  # Create fake data associated with that sample
  data %>%
    mutate(OD = y_rep_draw) -> draw_fake_data
  
  # Gets the parameters that generated this sample
  gather_draws(prior_samples$prior, cbind(gamma, kappa, epsilon)[k]) %>% 
    filter(.draw == draw) %>%
    add_column(distribution = "fake_data_draw") %>%
    mutate(Distribution = factor(distribution, levels = c("fake_data_draw"), labels = c("Sample from prior \nthat generates fake data"))) -> fake_data_check_params
  
  return(list(fake_data = draw_fake_data, true_params = fake_data_check_params))
}

#' get_prior_samples
#' 
#' This function takes a stan model and samples from the prior and returns all prior samples, the priors of the drug parameters and a list of 
#' samples that are appropriate for performing a fake data check with (non-flat).
#' 
#' @param stanmodel Stan object. Instance of a Stan model returned from stan_model() in Rstan.
#' @param stan_data List. List of data for Stan model.
#' @param no_of_chains Int. Number of chains for prior predictive sampling.
#' @param iter Int. Number of iterations for prior predictive sampling.
#' @param seed Int. Seed for sampling.
#' @param data Tibble. Tibble of real data.
#' 
#' @return List of prior samples, the priors for the drug parameters and the viable list of draws that can be used for fake data checks
#' 
get_prior_samples <- function(stanmodel, stan_data, no_of_chains, iter, seed, data) {
  # Sample from prior
  model_prior <- sampling(stanmodel, data = stan_data, chains = no_of_chains, iter = iter, seed = seed, cores = 1)
  
  # Get drug priors
  gather_draws(model_prior, cbind(gamma, kappa, epsilon)[k]) %>%
    add_column(distribution = "prior") -> priors
  
  # Locate divergent draws from priors (can happen with HS prior with no likelihood)
  prior_sampler_params <- get_sampler_params(model_prior, inc_warmup = FALSE)
  div_iters_list <- lapply(1:no_of_chains, function(x){tibble(.iteration = which(prior_sampler_params[[x]][, "divergent__"]==1), .chain=x)}) 
  do.call(rbind, div_iters_list) %>%
    add_column(divergent = TRUE) -> div_iters
  
  # Get draws from priors that are not divergent or flat. If a draw's trajectories are completely flat, drug-action parameters are non-identifiable by default.
  model_prior %>%
    spread_draws(y_tot[rowid]) %>%
    full_join(div_iters) %>%
    replace_na(list(divergent = FALSE)) %>%
    full_join(rowid_to_column(data), by = "rowid") %>%
    filter(time%in%c(min(data$time), max(data$time))) %>%
    group_by(.draw) %>%
    select(-c(rowid, .chain, .iteration, drug_conc, blanks, OD, drug1_name, drug2_name)) %>%
    pivot_wider(names_from = time, values_from = y_tot) %>%
    group_by(.draw, drug_conds, well) %>%
    mutate(growth_diff = `48.0053888888889` - `0`) %>%
    ungroup() %>%
    filter(drug_conds!="2X RPMI + dH2O") %>%  # filter out blanks
    select(.draw, drug_conds, growth_diff, divergent) %>%
    filter(!divergent) %>%
    group_by(.draw) %>%
    mutate(non_zero = growth_diff > 0) %>%  # only keep draws whose growth at end time point > initial inoculum
    summarise(keep = all(non_zero)) %>%
    filter(keep) %>%
    pull(.draw) -> draw_list
  
  return(list(prior = model_prior, drug_priors = priors, viable_draws = draw_list))
}

