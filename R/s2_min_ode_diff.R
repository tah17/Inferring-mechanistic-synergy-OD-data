rm(list = ls())
library(deSolve)
library(tidyverse)
library(tidybayes)
library(forecast)
library(cowplot)

# Functions --------------------------------------------------------------------
#
# Gompertz ODE with drug-action added
#
gompertz_drug <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX1 <- (beta/(1+gamma)) * X1 * log(K/(X1*(1+epsilon))) - kappa*X1 ## gomp with drug
    dX2 <- kappa*X1
    list(c(dX1, dX2))
  })
}
#
# Edwards ODE with drug-action added
#
edwards_drug <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX1 <- (beta/(1+gamma)) * X1 * exp(-(X1*(1+epsilon))/L) - kappa*X1 ## edwards with drug
    dX2 <- kappa*X1
    list(c(dX1, dX2))
  })
}
#
# Function that adds a small perturbation (perturb) to a chosen drug-action parameter (chosen_param), which can be either gamma, epsilon or kappa,
# and solves an ODE (ode_fn) at times, time, with parameters, params, and returns the ode solution (sol).
#
peturb_solve_ode <- function(ode_fn, time, params, perturb, chosen_param) {
  if (chosen_param == "gamma") {
    drug_params <- c(gamma=perturb, epsilon=0, kappa=0)
  } else if (chosen_param == "epsilon") {
    drug_params <- c(gamma=0, epsilon=perturb, kappa=0)
  } else {
    drug_params <- c(gamma=0, epsilon=0, kappa=perturb)
  }
  sol <- ode(y = c(X1 = params$y0, X2 = 0), times = time, func = ode_fn, parms = append(params, drug_params))
  return(sol)
}

# Solve ODEs --------------------------------------------------------------
#
# Read in ODE fits for Edwards and Gompertz models
#
gompertz_stan_fit <- readRDS(file=paste("output/gompertz_drug_hs_reg_full_data_full_fit.Rda", sep=""))  
edwards_stan_fit <- readRDS(file=paste("output/edwards_drug_hs_reg_full_data_full_fit.Rda", sep=""))  
#
# Get fitted Edwards and Gompertz models' parameter estimates at the mode
#
gather_draws(gompertz_stan_fit, cbind(beta, K, y0, basal, delta)) %>%
  group_by(.variable) %>%
  summarise(mode = Mode(.value)) %>%
  pivot_wider(names_from = .variable, values_from = mode) -> params_g

gather_draws(edwards_stan_fit, cbind(beta, L, y0, basal, delta)) %>%
  group_by(.variable) %>%
  summarise(mode = Mode(.value)) %>%
  pivot_wider(names_from = .variable, values_from = mode) -> params_e
#
#  Solve the ODEs when all drug-action parameters are at 0 (No drug solution)
#
time <- seq(0, 48, by = 1)
no_drug_edwards <- ode(y = c(X1 = params_e$y0, X2 = 0), times = time, func = edwards_drug, parms = append(params_e, c(gamma=0, epsilon=0, kappa=0)))
no_drug_gompertz <- ode(y = c(X1 = params_g$y0, X2 = 0), times = time, func = gompertz_drug, parms = append(params_g, c(gamma=0, epsilon=0, kappa=0)))
#
#  Solve both ODEs when each drug-action parameter is perturbed from 0 for a range of perturbations 
#
pertubations <- 10^(-6:-1)
parameters <- c("gamma", "epsilon", "kappa")  # drug parameters
pertubations_parameters <- expand.grid(pertubation = pertubations, parameter = parameters)
edwards_list <- lapply(1:nrow(pertubations_parameters), function(x) peturb_solve_ode(edwards_drug, time, params_e, pertubations_parameters$pertubation[x], pertubations_parameters$parameter[x]) %>% 
                        as.data.frame() %>%
                        tibble() %>% 
                        add_column(perturb=pertubations_parameters$pertubation[x]) %>% 
                        add_column(param=pertubations_parameters$parameter[x]) %>%
                        add_column(model="Edwards") %>%
                        mutate(y = params_e$basal + ((X1 + X2)/params_e$delta)) %>%  # get basal + (f_v + f_d)/delta (median model output)
                        mutate(RMSE = accuracy(params_e$basal + (no_drug_edwards[,"X1"] + no_drug_edwards[,"X2"])/params_e$delta, y)[, "RMSE"]) %>%  # RMSE between no drug and drug median model outputs
                        mutate(MAPE = accuracy(no_drug_edwards[,"X1"], X1)[, "MAPE"]))  # calculate percentage error between no drug and drug viable fungal growth 

gompertz_list <- lapply(1:nrow(pertubations_parameters), function(x) peturb_solve_ode(gompertz_drug, time, params_g, pertubations_parameters$pertubation[x], pertubations_parameters$parameter[x]) %>% 
                          as.data.frame() %>%
                          tibble() %>% 
                          add_column(perturb=pertubations_parameters$pertubation[x]) %>% 
                          add_column(param=pertubations_parameters$parameter[x]) %>%
                          add_column(model="Gompertz") %>%
                          mutate(y = params_g$basal + ((X1 + X2)/params_g$delta)) %>%  # get basal + (f_v + f_d)/delta (median model output)
                          mutate(RMSE = accuracy(params_g$basal + (no_drug_gompertz[,"X1"] + no_drug_gompertz[,"X2"])/params_g$delta, y)[, "RMSE"]) %>%  # RMSE between no drug and drug median model outputs
                          mutate(MAPE = accuracy(no_drug_gompertz[,"X1"], X1)[, "MAPE"]))  # calculate percentage error between no drug and drug viable fungal growth 
#
# combine and plot model solutions
#
do.call(rbind, edwards_list) %>%
  add_row(do.call(rbind, gompertz_list)) -> model_sols
#
# plot ODE solutions
#
model_sols %>%
  mutate(perturb = log10(perturb)) %>%
  mutate(`Perturbation size` = paste("10^", perturb, sep="")) %>%
  ggplot(aes(x = time, y = y)) +
  geom_line(aes(colour=param), size=0.5, alpha=0.9) +
  geom_line(data = tibble(model = "Edwards", time = no_drug_edwards[,"time"], y = params_e$basal + (no_drug_edwards[,"X1"]+no_drug_edwards[,"X2"])/params_e$delta), aes(x = time, y = y), size=1, colour="black", linetype="dashed") +
  geom_line(data = tibble(model = "Gompertz", time = no_drug_gompertz[,"time"], y = params_g$basal + (no_drug_gompertz[,"X1"]+no_drug_gompertz[,"X2"])/params_g$delta), aes(x = time, y = y), size=1, colour="black", linetype="dashed") +
  facet_grid(model~`Perturbation size`, scales="free", labeller="label_parsed") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.placement = "outside") +
  ylab("OD") +
  xlab("Time [hrs]") +
  scale_color_manual(name="Drug Parameter", breaks=c("gamma", "epsilon", "kappa"), values=c("#377eb8", "#4daf4a", "#e41a1c"), labels=c(expression(gamma), expression(epsilon), expression(kappa))) -> p_a
#
# plot MAPE and RMSE
#
model_sols %>%
  pivot_longer(RMSE:MAPE, names_to="metric") %>%
  select(perturb, param, model, metric, value) %>%
  unique() %>%
  ggplot(aes(x = perturb, y = value)) +
  geom_line(aes(colour=param), size=1) +
  facet_grid(metric ~ model, scales="free", switch="y") +
  geom_vline(xintercept = 1e-4, linetype="dashed", alpha=0.4) +
  geom_text(data = data.frame(metric = "RMSE", model="Gompertz", x = 1.6e-3, y = 5e-7, label = "Chosen threshold"), aes(x = x, y = y, label = label), size=3, colour="#636363") +
  geom_hline(data = data.frame(metric = "MAPE", yintercept = 0.1), aes(yintercept = yintercept), linetype="dashed") +
  geom_text(data = data.frame(metric = "MAPE", model="Edwards", x = 2e-6, y = 0.2, label = "0.1%"), aes(x = x, y = y, label = label), size=3.5) +
  geom_hline(data = data.frame(metric = "MAPE", yintercept = 1), aes(yintercept = yintercept), linetype="dashed") +
  geom_text(data = data.frame(metric = "MAPE", model="Edwards", x = 5e-6, y = 2, label = "1%"), aes(x = x, y = y, label = label), size=3.5) +
  geom_hline(data = data.frame(metric = "RMSE", yintercept = 0.009192358), aes(yintercept = yintercept), linetype="dashed") +
  geom_text(data = data.frame(metric = "RMSE", model="Edwards", x = 7.5e-6, y = 2e-2, label = "Lowest CV score"), aes(x = x, y = y, label = label), size=3) +
  scale_y_continuous(trans="log10", labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(trans="log10", breaks = pertubations, labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", strip.placement = "outside", strip.background = element_blank()) +
  xlab("Parameter perturbation") +
  ylab("") +
  scale_color_manual(name="Drug Parameter", breaks=c("gamma", "epsilon", "kappa"), values=c("#377eb8", "#4daf4a", "#e41a1c" ), labels=c(expression(gamma), expression(epsilon), expression(kappa))) -> p_b

legend <- get_legend(
  # create some space to the left of the legend
  p_b + theme(legend.box.margin = margin(0, 0, 0, 15), legend.key.size = unit(2, "line"), legend.position = "right") + guides(colour = guide_legend(override.aes = list(size=3.5)))
)

p_b_legend <- plot_grid(p_b + theme(legend.position = "none"), legend, rel_widths = c(0.72, 0.28), nrow = 1)
p <- plot_grid(p_a, p_b_legend, labels = "auto", rel_heights = c(0.8, 1), ncol = 1, label_size = 12) 

tiff("figures/figs4.tif", width = 20, height = 23, units = "cm", res=300)
p
dev.off()
