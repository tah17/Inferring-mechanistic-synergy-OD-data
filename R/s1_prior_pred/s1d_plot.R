rm(list = ls())
library(tidyverse)
library(tidybayes)
library(cowplot)
library(ggh4x)
seed <- 404806
set.seed(seed)

# Function ----------------------------------------------------------------
# Takes (1-alpha) level (desired_width), a tibble of parameters' posterior samples (posterior_params) and 
# the true parameter values (fake_data_params) and outputs the coverage for 
# each group of drug-action parameters: gamma, epsilon and kappa.
#
get_coverage_p_v_at_alpha <- function(desired_width, posterior_params, fake_data_params) {
  posterior_params %>%
    group_by(k, .variable, Draw, Model) %>%
    summarise(mode_hdi(.value, .width=desired_width)) %>%
    mutate(interval_id = row_number()) %>%   # some mode_hdi intervals can be discontinuous for some desired widths
    select(-c(y, .width, .point, .interval)) %>%
    full_join(fake_data_params) %>%
    group_by(k, .variable, Draw, Model, interval_id) %>%
    mutate(in_interval = .value >= ymin & .value <= ymax) %>%
    select(k, .variable, Draw, Model, interval_id, in_interval) %>%
    group_by(k, .variable, Draw, Model) %>%
    mutate(in_any_interval = any(in_interval)) %>%
    select(k, .variable, Draw, Model, in_any_interval) %>%
    unique() %>%
    group_by(Model, .variable) %>%
    summarise(raw_coverage = sum(in_any_interval), n = n(), coverage = sum(in_any_interval)/n()) %>%
    mutate(alpha = desired_width) -> coverage
  return(coverage)
}

# Read in results ---------------------------------------------------------
model_list <- c("edwards_drug_hs_reg_full_data", 
                "gompertz_drug_hs_reg_full_data")
n_draws <- 5
file_path <- "output/fake_data_check/"
draw_model_list <- expand.grid(Draw = 1:n_draws, Model = 1:length(model_list))  
# reads in fake data check results
list_priors <- lapply(model_list, function(x){readRDS(file=paste(file_path, "drug_priors_", x, ".Rda", sep = "")) %>% add_column(Model = x)})
priors <- do.call(rbind, list_priors)
list_fake_data_params <- lapply(1:nrow(draw_model_list), function(x){readRDS(file=paste(file_path, "fake_data_params_", x, ".Rda", sep = "")) %>% add_column(Model = model_list[draw_model_list$Model[x]]) %>% add_column(Draw = draw_model_list$Draw[x])})
fake_data_params <- do.call(rbind, list_fake_data_params)
list_posterior_params <- lapply(1:nrow(draw_model_list), function(x){readRDS(file=paste(file_path, model_list[draw_model_list$Model[x]], "_posterior_params_", draw_model_list$Draw[x], ".Rda", sep = "")) %>% add_column(distribution = "posterior") %>% add_column(Model = model_list[draw_model_list$Model[x]]) %>% add_column(Draw = draw_model_list$Draw[x])})
posterior_params <- do.call(rbind, list_posterior_params)

lapply(1:n_draws, function(x){add_column(priors, Draw=x)}) %>%
  do.call(rbind, .) -> prior_params

rbind(prior_params, posterior_params, fake_data_params) %>%
  mutate(Distribution = factor(distribution, levels = c("prior", "posterior", "Sample from prior \nthat generates fake data"), labels = c("Prior", "Posterior \n(fake data check)", "Sample from prior \nthat generates fake data"))) %>%
  select(-c(.chain, .iteration, .draw)) -> draws

# Calculate posterior contraction -----------------------------------------
#
#  gets z-scores of fake data check
#
posterior_params %>%
  group_by(k, .variable, Draw, Model) %>%
  reframe(mu_post = mean(.value), sd_post = sd(.value)) %>%
  full_join(fake_data_params) %>%
  group_by(k, .variable, Draw, Model) %>%
  mutate(z = (mu_post - .value)/sd_post) %>%
  select(k, .variable, Draw, Model, z, sd_post) -> z_scores
#
#  gets s-scores of fake data check
#
prior_params %>%
  group_by(k, .variable, Draw, Model) %>%
  reframe(sd_prior = sd(.value)) %>%
  full_join(z_scores) %>%
  group_by(k, .variable, Draw, Model) %>%
  mutate(s = 1-(sd_post^2)/(sd_prior^2)) %>%
  select(k, .variable, Draw, Model, z, s) -> z_scores_tot
#
#  plots s, z-scores
#
z_scores_tot %>%
  mutate(Model = factor(Model, levels = model_list, labels = c("Edwards-D-HS", "Gompertz-D-HS"))) %>%
  mutate(k = case_when(k==1 ~ "AFG", k==2 ~ "MGX", k==3 ~ "5FC", k==4 ~ "AFG:MGX", k==5 ~ "AFG:5FC")) %>%
  mutate(k = factor(k, levels = c("AFG", "MGX", "5FC", "AFG:MGX", "AFG:5FC"))) %>%
  mutate(action = case_when(k%in%c("AFG", "MGX", "5FC") ~ "Individual", k%in%c("AFG:MGX", "AFG:5FC") ~ "Synergy")) %>%
  mutate(.variable = factor(.variable, levels = c("gamma", "epsilon", "kappa"))) %>%
  mutate(k_x_axis = paste(.variable, k, sep="_")) %>%
  mutate(k_x_axis = factor(k_x_axis, levels=c("gamma_AFG", "gamma_MGX", "gamma_5FC", "gamma_AFG:MGX", "gamma_AFG:5FC", "epsilon_AFG", "epsilon_MGX", "epsilon_5FC", "epsilon_AFG:MGX", "epsilon_AFG:5FC", "kappa_AFG", "kappa_MGX", "kappa_5FC", "kappa_AFG:MGX", "kappa_AFG:5FC"))) %>%
  ggplot(aes(x = s, y = z, color=.variable)) +
  geom_vline(xintercept = 1, linetype="dashed", alpha=0.8) +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.8) +
  geom_hline(yintercept = 3, linetype="dashed", alpha=0.4) +
  geom_hline(yintercept = -3, linetype="dashed", alpha=0.4) +
  geom_point(size = 2, stroke=1, shape=22) +
  facet_grid(Model ~ Draw) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", axis.text.x = element_text(angle=25, vjust=1, hjust=1)) +
  xlim(0, 1) +
  ylim(-4, 4) +
  scale_color_manual(name="Drug Parameter", breaks=c("gamma", "epsilon", "kappa"), values=c("#377eb8", "#4daf4a", "#e41a1c" ), labels=c(expression(gamma), expression(epsilon), expression(kappa))) -> p_a
#
# gets distance to ideal (s, z)-score of (1, 0)
#
z_scores_tot %>%
  group_by(k, .variable, Draw, Model) %>%
  summarise(dist_to_ideal = dist(rbind(c(1,0), c(s, z)))) -> dists
#
# plots distances
#
model_colours <- c("#4d004b", "#045a8d")
dists %>%
  ggplot(aes(x = dist_to_ideal, fill=Model, color=Model)) +
  geom_density(adjust=1.5, alpha=0.3) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  ylab("Density") +
  xlab("Distance of (S, Z)-scores to Ideal") +
  scale_fill_manual(name="Model", breaks=model_list, values=model_colours, labels=c("Edwards-D-HS", "Gompertz-D-HS")) +
  scale_color_manual(name="Model", breaks=model_list, values=model_colours, labels=c("Edwards-D-HS", "Gompertz-D-HS")) -> p_b
#
#  calculates and plots coverage 
#
tot_coverage_p_v <- lapply(seq(0.1, 1, length.out=10), function(x){get_coverage_p_v_at_alpha(x, posterior_params, fake_data_params)})
do.call(rbind, tot_coverage_p_v) %>%
  mutate(Model = factor(Model, levels = model_list, labels = c("Edwards-D-HS", "Gompertz-D-HS"))) %>%
  rowwise() %>%
  mutate(binom_ci_l = prop.test(raw_coverage, n)$conf.int[1], binom_ci_u = prop.test(raw_coverage, n)$conf.int[2]) %>%
  ggplot(aes(x = alpha, y = coverage, group=interaction(Model, .variable), fill=.variable)) +
  facet_grid(. ~ Model) +
  geom_ribbon(aes(ymin = binom_ci_l, ymax = binom_ci_u), alpha = 0.1) +
  geom_line(aes(color=.variable), size=1) +
  geom_abline(intercept=0, slope = 1, linetype = "dashed", size=1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top") +
  ylab("Coverage") +
  xlab("1 - alpha") +
  scale_color_manual(name="Drug Parameter", breaks=c("gamma", "epsilon", "kappa"), values=c("#377eb8", "#4daf4a", "#e41a1c" ), labels=c(expression(gamma), expression(epsilon), expression(kappa))) +
  scale_fill_manual(name="Drug Parameter", breaks=c("gamma", "epsilon", "kappa"), values=c("#377eb8", "#4daf4a", "#e41a1c" ), labels=c(expression(gamma), expression(epsilon), expression(kappa))) -> p_c

p1 <- plot_grid(p_b + theme(legend.position = "top"), p_c + theme(legend.position = "top"), labels = c('b', 'c'), rel_widths = c(0.8, 1), nrow = 1, label_size = 12) 
p2 <- plot_grid(p_a, p1, labels = c('a', ''), rel_heights = c(1, 0.9), ncol = 1, label_size = 12) 
tiff("figures/figs2.tif", width = 19, height = 17, units = "cm", res=300)
p2
dev.off()

# Plot Fake Data Check ----------------------------------------------------
chosen_model <- 1  # the edwards_drug_hs_reg_full_data model
rbind(prior_params, posterior_params, fake_data_params) %>%
  mutate(Distribution = factor(distribution, levels = c("prior", "posterior", "fake_data_draw"), labels = c("Prior", "Posterior \n(fake data check)", "Sample from prior \nthat generates fake data"))) %>%
  select(-c(.chain, .iteration, .draw)) %>%
  mutate(k = case_when(k==1 ~ "AFG", k==2 ~ "MGX", k==3 ~ "5FC", k==4 ~ "AFG:MGX", k==5 ~ "AFG:5FC")) %>%
  mutate(k = factor(k, levels = c("AFG", "MGX", "5FC", "AFG:MGX", "AFG:5FC"))) %>%
  mutate(.variable = factor(.variable, levels = c("gamma", "epsilon", "kappa"))) -> draws
ggplot() +
  stat_pointinterval(data = filter(draws, Model==model_list[chosen_model], !distribution%in%c("fake_data_draw")), aes(x = k, y = .value, color=Distribution), position="dodge", point_interval = mode_hdi, .width = c(0.80, 0.95), point_size=1.5) +
  geom_point(data = filter(draws, Model==model_list[chosen_model], distribution=="fake_data_draw"), aes(y = .value, x = k, shape=Distribution), size=1.5, stroke=1, fill="white") +
  theme_bw(base_size = 11) +
  scale_color_manual(values=c("#cab2d6", "#6a3d9a")) +
  facet_grid(.variable ~ Draw, scales="free", labeller="label_parsed") +
  scale_shape_manual(name = "", values=21) +
  xlab("Drug condition") +
  ylab("Value") +
  theme(legend.position="top", axis.text.x = element_text(angle=25, vjust=1, hjust=1)) -> p3

tiff("figures/figs3.tif", width = 19, height = 13, units = "cm", res=300)
p3
dev.off()

# to calculate the coverage of the parameter estimates shown in p3 run the below:
# get_coverage_p_v_at_alpha(0.95, posterior_params, fake_data_params)
  