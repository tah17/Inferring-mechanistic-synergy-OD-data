rm(list = ls())
library(ggplot2)
library(ggdist)
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(rstan)
library(cowplot)
library(ggh4x)
seed <- 404806
set.seed(seed)

# Read in Fitted Model ----------------------------------------------------
model <- "gompertz_drug_hs_reg_full_data"
stan_fit <- readRDS(file=paste("output/", model, "_full_fit.Rda", sep=""))  #reads model
#
#  plot inferred drug action parameters
#
pal_fill <- scales::brewer_pal(palette = "Dark2")(5)  # colour palette
label_map <- c(
  "gamma_AFG" = expression(gamma[AFG]),
  "gamma_MGX" = expression(gamma[MGX]),
  "gamma_5FC" = expression(gamma[`5FC`]),
  "gamma_AFG:MGX" = expression(gamma[AFG:MGX]),
  "gamma_AFG:5FC" = expression(gamma[AFG:`5FC`]),
  "epsilon_AFG" = expression(epsilon[AFG]),
  "epsilon_MGX" = expression(epsilon[MGX]),
  "epsilon_5FC" = expression(epsilon[`5FC`]),
  "epsilon_AFG:MGX" = expression(epsilon[AFG:MGX]),
  "epsilon_AFG:5FC" = expression(epsilon[AFG:`5FC`]),
  "kappa_AFG" = expression(kappa[AFG]),
  "kappa_MGX" = expression(kappa[MGX]),
  "kappa_5FC" = expression(kappa[`5FC`]),
  "kappa_AFG:MGX" = expression(kappa[AFG:MGX]),
  "kappa_AFG:5FC" = expression(kappa[AFG:`5FC`])
)
gather_draws(stan_fit, cbind(gamma, kappa, epsilon)[k]) %>%
  mutate(k = case_when(k==1 ~ "AFG", k==2 ~ "MGX", k==3 ~ "5FC", k==4 ~ "AFG:MGX", k==5 ~ "AFG:5FC")) %>%
  mutate(k = factor(k, levels = c("AFG", "MGX", "5FC", "AFG:MGX", "AFG:5FC"))) %>%
  mutate(action = case_when(k%in%c("AFG", "MGX", "5FC") ~ "Individual", k%in%c("AFG:MGX", "AFG:5FC") ~ "Synergy")) %>%
  mutate(.variable = factor(.variable, levels = c("gamma", "epsilon", "kappa"))) %>%
  mutate(k_x_axis = paste(.variable, k, sep="_")) %>%
  mutate(k_x_axis = factor(k_x_axis, levels=c("gamma_AFG", "gamma_MGX", "gamma_5FC", "gamma_AFG:MGX", "gamma_AFG:5FC", "epsilon_AFG", "epsilon_MGX", "epsilon_5FC", "epsilon_AFG:MGX", "epsilon_AFG:5FC", "kappa_AFG", "kappa_MGX", "kappa_5FC", "kappa_AFG:MGX", "kappa_AFG:5FC"))) %>%
  mutate(variable_fn = factor(.variable, levels = c("gamma", "epsilon", "kappa"), labels=c("Inhibition of Growth Rate", "Reduction of Carrying Capacity",  "Direct Killing"))) %>%
  group_by(k_x_axis) %>%
  mutate(ymin = min(mode_hdi(.value, .width=0.95)$ymin)) %>%  
  mutate(non_zero = case_when(ymin > 1e-4 ~ TRUE, .default = FALSE)) %>%  # only shade in parameters whose 95% CIs lower bound > 10^(-4)
  ungroup() %>%
  mutate(colour_for_plot = case_when(!non_zero ~ "NA", non_zero&k=="AFG" ~ "AFG", non_zero&k=="MGX" ~ "MGX", non_zero&k=="5FC" ~ "5FC", non_zero&k=="AFG:MGX" ~ "AFG+MGX", non_zero&k=="AFG:5FC" ~ "AFG+5FC")) %>%
  ggplot() +
  stat_pointinterval(aes(x=k_x_axis, y = .value, colour = colour_for_plot), point_interval = mode_hdi, .width = c(0.8, 0.95)) +
  theme_bw(base_size = 11) +
  facet_nested_wrap(~ variable_fn + action, scales="free", nrow=1) +
  scale_x_discrete(labels = label_map) +
  theme(legend.position="top", axis.text.x = element_text(angle=25, vjust=1, hjust=1)) +
  xlab("Drug parameters") +
  ylab("Parameter value") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(name="Acting Drug Condition", breaks=c("NA", "AFG", "MGX", "5FC", "AFG+MGX", "AFG+5FC"), values=c("#969696", pal_fill)) -> p

tiff("figures/figs5.tif", width = 20, height = 13, units = "cm", res=300)
p
dev.off()
