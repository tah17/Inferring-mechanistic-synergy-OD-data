rm(list = ls())
library(tidyverse)
library(tidybayes)
library(cowplot)
library(ggh4x)
seed <- 404806
set.seed(seed)

# Plot Fake Data Check ----------------------------------------------------
model <- "edwards_drug_hs_reg_fixed_k_full_data"
fake_data_params <- readRDS(file=paste("output/fake_data_check/high_kappa/fake_data_params_", model, ".Rda", sep=""))  
posterior_params <- readRDS(paste("output/fake_data_check/high_kappa/", model, "_posterior_params", ".Rda", sep = "")) %>% add_column(distribution="posterior")
prior <- readRDS(file=paste("output/fake_data_check/drug_priors_edwards_drug_hs_reg_full_data.Rda", sep = "")) 

rbind(prior, posterior_params, fake_data_params) %>%
  mutate(Distribution = factor(distribution, levels = c("prior", "posterior", "fake_data_draw"), labels = c("Prior", "Posterior \n(fake data check)", "High Kappa Sample \nthat generates fake data"))) %>%
  select(-c(.chain, .iteration, .draw)) %>%
  mutate(k = case_when(k==1 ~ "AFG", k==2 ~ "MGX", k==3 ~ "5FC", k==4 ~ "AFG:MGX", k==5 ~ "AFG:5FC")) %>%
  mutate(k = factor(k, levels = c("AFG", "MGX", "5FC", "AFG:MGX", "AFG:5FC"))) %>%
  mutate(.variable = factor(.variable, levels = c("gamma", "epsilon", "kappa"))) -> draws
ggplot() +
  stat_pointinterval(data = filter(draws, !distribution%in%c("fake_data_draw")), aes(x = k, y = .value, color=Distribution), position="dodge", point_interval = mode_hdi, .width = c(0.80, 0.95), point_size=1.5) +
  geom_point(data = filter(draws, distribution=="fake_data_draw"), aes(y = .value, x = k, shape=Distribution), size=1.5, stroke=1, fill="white") +
  theme_bw(base_size = 11) +
  scale_color_manual(values=c("#cab2d6", "#6a3d9a")) +
  facet_wrap(. ~ .variable, scales="free", labeller="label_parsed") +
  scale_shape_manual(name = "", values=21) +
  xlab("Drug condition") +
  ylab("Value") +
  theme(legend.position="top", axis.text.x = element_text(angle=25, vjust=1, hjust=1)) -> p

tiff("figures/figS5.tif", width = 19, height = 10, units = "cm", res=300)
p
dev.off()
