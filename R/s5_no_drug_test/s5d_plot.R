rm(list = ls())
library(ggplot2)
library(ggdist)
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(rstan)
library(cowplot)
seed <- 404806
set.seed(seed)

# Read in Fitted Model ----------------------------------------------------
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
no_of_trials <- 3
file_path <- "output/no_drug_test/"
draw_model_list <- expand.grid(Trial = 1:no_of_trials, Model = 1:length(model_list))  
# reads in fake data check results
list_posterior_params <- lapply(1:nrow(draw_model_list), function(x){readRDS(file=paste(file_path, model_list[draw_model_list$Model[x]], "_posterior_params_", draw_model_list$Trial[x], ".Rda", sep = "")) %>% add_column(Model = model_list[draw_model_list$Model[x]]) %>% add_column(Trial = draw_model_list$Trial[x])})
posterior_params <- do.call(rbind, list_posterior_params)
#
#  plot inferred drug action parameters
#
label_map <- c(
  "gamma_test" = expression(gamma[test]),
  "epsilon_test" = expression(epsilon[test]),
  "kappa_test" = expression(kappa[test])
)
posterior_params %>%
  mutate(Model = factor(Model, levels = model_list, labels=c("Edwards-D-HS", "Gompertz-D-HS"))) %>%
  mutate(Trial = factor(Trial, levels = 1:3, labels=paste("Rep", 1:3, sep=" "))) %>%
  mutate(k = case_when(k==1 ~ "test")) %>%
  mutate(k = factor(k)) %>%
  mutate(.variable = factor(.variable, levels = c("gamma", "epsilon", "kappa"))) %>%
  mutate(k_x_axis = paste(.variable, k, sep="_")) %>%
  mutate(k_x_axis = factor(k_x_axis, levels=c("gamma_test","epsilon_test", "kappa_test"))) %>%
  group_by(k_x_axis, Model, Trial) %>%
  mutate(ymin = min(mode_hdi(.value, .width=0.95)$ymin)) %>%
  mutate(non_zero = case_when(ymin > 1e-4 ~ TRUE, .default = FALSE)) %>%  # only shade in parameters whose 95% CIs lower bound > 10^(-4)
  ungroup() %>%
  mutate(colour_for_plot = case_when(!non_zero ~ "NA", non_zero&k=="test" ~ "test")) %>%
  ggplot() +
  stat_pointinterval(aes(x=k_x_axis, y = .value, colour = colour_for_plot), point_interval = mode_hdi, .width = c(0.8, 0.95)) +
  theme_bw(base_size = 11) +
  facet_grid(Trial ~ Model) +
  scale_x_discrete(labels = label_map) +
  theme(legend.position="none") +
  xlab("Drug parameters") +
  ylab("Parameter value") +
  scale_color_manual(breaks=c("NA", "test"), values=c("#969696", "red")) -> p

tiff("figures/figS12.tif", width = 20, height = 20, units = "cm", res=300)
p
dev.off()

