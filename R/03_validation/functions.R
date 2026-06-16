library(forecast)
library(loo)
library(tidybayes)
#' get_metrics
#' 
#' This function takes a stan fit and data labelled as training or testing and returns fit + prediction metrics
#' 
#' @param dat_i Tibble. Data labelled as training or testing of fold i..
#' @param stan_i Stan object. Stan fit object trained on training data of fold i.
#' 
#' @return Tibble of Mean and SEM of ME, RMSE, MAE, MPE, MAPE & LPD
#' 
get_metrics <- function(dat_i, stan_i) {
  # calculates forecast::accuracy between model fit/predictions and  train/test data (resp.)
  stan_i %>%
    spread_draws(y_tot[rowid]) %>%
    full_join(rowid_to_column(dat_i), by = "rowid") %>%
    group_by(well, .draw, testing) %>%  # calculates metrics per draw and well
    reframe(accs = as.list(accuracy(OD, y_tot))) %>%
    group_by(well, .draw, testing) %>%  # reframe returns ungrouped tibble
    mutate(accs = as.numeric(accs)) %>%
    mutate(metric = c('ME', 'RMSE', 'MAE', 'MPE', 'MAPE')) %>%  # metrics in forecast::accuracy
    group_by(well, metric, testing) %>%  
    mutate(accs = mean(accs)) %>%  # get mean metrics over posterior draws per well
    ungroup() %>%
    select(-.draw) %>%
    unique() %>%
    group_by(metric, testing) %>%  # averages metrics over replicates in the wells
    summarise(mean = mean(accs), sem = sd(accs)/sqrt(n()), .groups="keep") %>%
    distinct(metric, testing, mean, sem) %>%
    ungroup() -> accs_metrics

  log_lik <- extract_log_lik(stan_i)  # extracts log likelihoods
  lpds <- colMeans(exp(log_lik))  # calcs average likelihood over Stan iters 
  # calculates lpds at test data 
  dat_i %>%
    filter(testing) %>%
    mutate(lpd = log(lpds)) %>%  # converts to log likelihoods
    group_by(well) %>%  # calculates mean lpd per replicate in a well
    summarise(mean_lpd = mean(lpd)) %>%
    mutate(mean_lpd = as.numeric(mean_lpd)) %>%
    ungroup() %>%
    summarise(mean = mean(mean_lpd), sem = sd(mean_lpd)/sqrt(n())) %>%  # averages lpds over replicates in the wells
    distinct(mean, sem) %>%
    ungroup() %>%
    add_column(metric = "LPD") %>%
    add_column(testing = TRUE) -> lpd_metrics
  
  # combines forecast accuracies and lpds
  accs_metrics %>%
    add_row(lpd_metrics) -> metrics

  return(metrics)
}