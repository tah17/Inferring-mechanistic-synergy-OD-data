rm(list = ls())
library(ggplot2)
library(ggdist)
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(rstan)
library(cowplot)
library(ggh4x)
library(GGally)
seed <- 404806
set.seed(seed)

# Read in Fitted Model ----------------------------------------------------
model <- "edwards_drug_hs_reg_full_data"
stan_fit <- readRDS(file=paste("output/", model, "_full_fit.Rda", sep=""))  #reads model
#
# scatter plot of inferred drug action parameters
#
label_map <- c(
  "gamma_AFG" = "gamma[AFG]",
  "gamma_MGX" = "gamma[MGX]",
  "gamma_5FC" = "gamma[`5FC`]",
  "gamma_AFG:MGX" = "gamma[AFG:MGX]",
  "gamma_AFG:5FC" = "gamma[AFG:`5FC`]",
  "epsilon_AFG" = "epsilon[AFG]",
  "epsilon_MGX" = "epsilon[MGX]",
  "epsilon_5FC" = "epsilon[`5FC`]",
  "epsilon_AFG:MGX" = "epsilon[AFG:MGX]",
  "epsilon_AFG:5FC" = "epsilon[AFG:`5FC`]",
  "kappa_AFG" = "kappa[AFG]",
  "kappa_MGX" = "kappa[MGX]",
  "kappa_5FC" = "kappa[`5FC`]",
  "kappa_AFG:MGX" = "kappa[AFG:MGX]",
  "kappa_AFG:5FC" = "kappa[AFG:`5FC`]"
)
gather_draws(stan_fit, cbind(gamma, kappa, epsilon)[k]) %>%
  mutate(k = case_when(k==1 ~ "AFG", k==2 ~ "MGX", k==3 ~ "5FC", k==4 ~ "AFG:MGX", k==5 ~ "AFG:5FC")) %>%
  mutate(k = factor(k, levels = c("AFG", "MGX", "5FC", "AFG:MGX", "AFG:5FC"))) %>%
  mutate(.variable = factor(.variable, levels = c("gamma", "epsilon", "kappa"))) %>%
  mutate(k_x_axis = paste(.variable, k, sep="_")) %>%
  mutate(k_x_axis = factor(k_x_axis, levels=c("gamma_AFG", "gamma_MGX", "gamma_5FC", "gamma_AFG:MGX", "gamma_AFG:5FC", "epsilon_AFG", "epsilon_MGX", "epsilon_5FC", "epsilon_AFG:MGX", "epsilon_AFG:5FC", "kappa_AFG", "kappa_MGX", "kappa_5FC", "kappa_AFG:MGX", "kappa_AFG:5FC"))) %>%
  ungroup() %>%
  select(-c(.chain, .iteration, k, .variable)) %>%
  pivot_wider(names_from = k_x_axis, values_from = .value) %>%
  ggpairs(columns = 2:16, 
          labeller = as_labeller(label_map, default = label_parsed),
          lower = list(continuous = wrap("density", colour = "grey40", fill = "grey80"))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) -> p

tiff("figures/figS7.tif", width = 30, height = 30, units = "cm", res=300)
p
dev.off()
  