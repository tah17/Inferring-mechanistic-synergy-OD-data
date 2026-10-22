#
# Script that uses posterior predictive to visualise different MGX application times 
#
rm(list = ls())
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(deSolve)
seed <- 404806
set.seed(seed)

# Edwards ODE -------------------------------------------------------------
#
# Edwards ODE with drug-action added
#
edwards_drug <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX1 <- (beta/(1+gamma)) * X1 * exp(-(X1*(1+epsilon))/L) - kappa*X1  # edwards with drug
    dX2 <- kappa*X1
    list(c(dX1, dX2))
  })
}

# Specify and read in model ---------------------------------------------------
model <- "edwards_drug_hs_reg_full_data"
fit <- readRDS(file=paste("output/", model, "_full_fit.Rda", sep=""))  #reads model

# Get Posterior Parameters ------------------------------------------------
#
# Get median parameter values over draws
#
gather_draws(fit, cbind(gamma, kappa, epsilon)[k]) %>%
  unite(.variable, c(.variable, k)) %>%
  full_join(gather_draws(fit, cbind(beta, L, y0))) %>%
  group_by(.variable) %>%
  summarise(median = median(.value)) %>%
  pivot_wider(names_from = .variable, values_from = median) -> post_params
#
# Get AFG drug-action median parameters to act as a baseline before MGX administration
#
post_params %>%
  select(L, beta, epsilon_1, gamma_1, kappa_1, y0) %>%
  rename(epsilon = epsilon_1, gamma = gamma_1, kappa = kappa_1) -> base_params
#
# AFG+MGX drug-action median parameters 
#
post_params %>%
  select(L, beta, epsilon_1, gamma_1, kappa_1, epsilon_2, gamma_2, kappa_2, epsilon_4, gamma_4, kappa_4, y0) %>%
  mutate(epsilon = epsilon_1 + epsilon_2 + epsilon_4, gamma = gamma_1 + gamma_2 + gamma_4, kappa = kappa_1 + kappa_2 + kappa_4) %>%
  select(L, beta, y0, epsilon, gamma, kappa) -> combo_drug_params

# Solve ODE at median params ----------------------------------------------
time <- seq(0, 48, by = 1) # time span of experiment
#
# No drug solution 
#
edwards_no_drug <- ode(y = c(X1 = base_params$y0, X2 = 0), times = time, func = edwards_drug, parms = c(beta = base_params$beta, L = base_params$L, epsilon = 0, kappa = 0, gamma = 0)) %>% 
  as_tibble() %>%
  add_column(drug_cond = "RPMI")
#
# AFG solution 
#
edwards_afg <- ode(y = c(X1 = base_params$y0, X2 = 0), times = time, func = edwards_drug, parms = base_params) %>%
  as_tibble() %>%
  add_column(drug_cond = "AFG")
#
# Combine both solutions
#
rbind(edwards_no_drug, edwards_afg) %>%
  mutate(X1 = as.numeric(X1)) %>%
  mutate(X2 = as.numeric(X2)) %>%
  mutate(time = as.numeric(time)) -> med_preds
#
# Get max value of f_v growing in RPMI at the end of the experiment for comparator used in plot
#
max_f <- med_preds %>%
  filter(time == max(time), drug_cond=="RPMI") %>%
  pull(X1)

# Solve ODE using 100 samples from posterior ------------------------------
#
# Sample 100 parameters from full posterior
#
gather_draws(fit, cbind(gamma, kappa, epsilon)[k]) %>%
  unite(.variable, c(.variable, k)) %>%
  full_join(gather_draws(fit, cbind(beta, L, y0))) %>%
  group_by(.variable) %>%
  pivot_wider(names_from = .variable, values_from = .value) -> full_params
sample_draws <- sample(1:nrow(full_params), 100)
full_params %>%
  filter(.draw %in% sample_draws) -> sampled_params
#
# Solve ODE for AFG sole administration using 100 samples from posterior
#
drug_preds <- lapply(1:nrow(sampled_params), function(x) ode(y = c(X1 = sampled_params$y0[x], X2 = 0), times = time, func = edwards_drug, parms = c(beta = sampled_params$beta[x], 
                                                                                                                                                    L = sampled_params$L[x], 
                                                                                                                                                    epsilon = sampled_params$epsilon_1[x], 
                                                                                                                                                    kappa = sampled_params$kappa_1[x], 
                                                                                                                                                    gamma = sampled_params$gamma_1[x])) %>%
                       as_tibble() %>%
                       add_column(drug_cond = "AFG") %>%
                       add_column(draw = x))

do.call(rbind, drug_preds) %>%
  mutate(X1 = as.numeric(X1)) %>%
  mutate(X2 = as.numeric(X2)) %>%
  mutate(time = as.numeric(time)) -> preds

# Solve ODE when adding diff drug times ---------------------------
added_times <- seq(from = 1, to = 45, by = 8)  # different times MGX is added
#
# Solve ODE at median parameter values when AFG is added at time 0 and MGX is added at different times 
#
med_added_drug_preds <- lapply(added_times, function(x) ode(y = c(X1 = filter(med_preds, drug_cond=="AFG", time == x)$X1, X2 = filter(med_preds, drug_cond=="AFG", time == x)$X2), times = time[which(time >= x)], func = edwards_drug, parms = combo_drug_params) %>%
                                 as_tibble() %>%
                                 rbind(edwards_afg[which(time < x) , ] %>% select(time:X2), .) %>%
                                 mutate(added_time = x))
med_added_preds <- do.call(rbind, med_added_drug_preds) %>% add_column(drug_cond = "AFG+MGX") 
#
# Check when median ODE of AFG + MGX is 90% below max f_v value at the end of the experiment
#
med_added_preds %>%
  group_by(added_time) %>%
  mutate(reduction_below_90 = case_when(time == max(time) & X1 < 0.1*max_f ~ TRUE, .default = FALSE)) %>%
  mutate(inhibited = any(reduction_below_90)) %>%
  ungroup() -> med_added_preds_labelled
#
# Solve ODE using 100 samples from posterior when AFG is added at time 0 and MGX is added at different times 
#
samples_times <- expand.grid(params = 1:nrow(sampled_params), times = added_times)
added_drug_preds <- lapply(1:nrow(samples_times), function(x) ode(y = c(X1 = filter(preds, draw == samples_times$params[x], drug_cond == "AFG", time == samples_times$times[x])$X1, X2 = filter(preds, draw == samples_times$params[x], drug_cond == "AFG", time == samples_times$times[x])$X2), times = time[which(time >= samples_times$times[x])], func = edwards_drug, parms = c(beta = sampled_params$beta[samples_times$params[x]], 
                                                                                                                                                          L = sampled_params$L[samples_times$params[x]], 
                                                                                                                                                          epsilon = sampled_params$epsilon_1[samples_times$params[x]] + sampled_params$epsilon_2[samples_times$params[x]] + sampled_params$epsilon_4[samples_times$params[x]], 
                                                                                                                                                          kappa = sampled_params$kappa_1[samples_times$params[x]] + sampled_params$kappa_2[samples_times$params[x]] + sampled_params$kappa_4[samples_times$params[x]], 
                                                                                                                                                          gamma = sampled_params$gamma_1[samples_times$params[x]] + sampled_params$gamma_2[samples_times$params[x]] + sampled_params$gamma_4[samples_times$params[x]])) %>%
                             as_tibble() %>%
                             rbind(filter(preds, draw == samples_times$params[x], drug_cond == "AFG", time < samples_times$times[x]) %>% select(time:X2), .) %>%
                             mutate(added_time = samples_times$times[x]) %>%
                             add_column(draw = samples_times$params[x]))
added_preds <- do.call(rbind, added_drug_preds) %>% add_column(drug_cond = "AFG+MGX") 

# Plot different adding drug times ----------------------------------------
added_preds %>%
  group_by(added_time, draw) %>%
  mutate(reduction_below_90 = case_when(time == max(time) & X1 < 0.1*max_f ~ TRUE, .default = FALSE)) %>%
  mutate(inhibited = any(reduction_below_90)) %>%
  ungroup() %>%
  mutate(drug_cond = case_when(time < added_time ~ "AFG", .default = drug_cond)) %>%
  ggplot(aes(x = time)) +
    geom_rect(data = tibble(added_time = added_times), aes(xmin=added_time, xmax=Inf, ymin=0, ymax=Inf), inherit.aes = FALSE, fill = "#e7298a", alpha = 0.1) +
    geom_rect(data = tibble(added_time = added_times), aes(xmin=0, xmax=added_time, ymin=0, ymax=Inf), inherit.aes = FALSE, fill = "#1b9e77", alpha = 0.1) +
    geom_line(aes(group=draw, y = X1), alpha = 0.1, linewidth = 0.5) +
    geom_line(data = med_added_preds_labelled, aes(y = X1, color = inhibited), linewidth = 1) +
    geom_vline(aes(xintercept = added_time), alpha = 0.8, linewidth = 0.75, linetype="dashed") +
    facet_wrap(. ~ added_time) +
    theme_bw(base_size = 11) +
    ylab(expression(f[v]~"[cells/"*mu*"l]")) +
    xlab("Time [hrs]") +
    theme(legend.position="top") +
    scale_y_continuous(trans="log10", labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(name="90% Growth Reduction", breaks = c(TRUE, FALSE), values=c("#1a1a1a", "#ae017e")) -> p

tiff("figures/figS11.tif", width = 15, height = 10, units = "cm", res=300)
p
dev.off()

