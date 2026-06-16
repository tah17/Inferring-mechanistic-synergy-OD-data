rm(list = ls())
library(tidyverse)
library(cowplot)
library(ggh4x)
seed <- 404806
set.seed(seed)
#
# specify all models used in CV: filenames, model names and model colours in plots
#
drug_free_model_list <- c("exponential", 
                          "logistic",
                          "gompertz",
                          "edwards",
                          "gompertz_OD_direct",
                          "edwards_OD_direct")
indv_drug_model_list <- c("gompertz_indv_drug",
                          "edwards_indv_drug", 
                          "edwards_indv_drug_gamma0", 
                          "edwards_indv_drug_epsilon0",
                          "edwards_indv_drug_kappa0",
                          "edwards_indv_drug_decay",
                          "edwards_indv_drug_hs_reg")
all_drug_model_list <- c("edwards_drug_hs_reg", 
                         "edwards_drug_no_syn_hs_reg",
                         "gompertz_drug_hs_reg",
                         "gompertz_drug_no_syn_hs_reg")
model_lists <- list(drug_free_model_list, indv_drug_model_list, all_drug_model_list)
conditions <- c("drug_free", "indv_drug", "all_drug")

model_names <- c("Exponential", 
                 "Logistic",
                 "Gompertz", 
                 "Edwards", 
                 "Gompertz\n(OD direct)",
                 "Edwards\n(OD direct)", 
                 "Gompertz-invD", 
                 "Edwards-indvD", 
                 "Edwards-indvD\n(Gamma=0)",
                 "Edwards-indvD\n(Epsilon=0)",
                 "Edwards-indvD\n(Kappa=0)",
                 "Edwards-indvD\n(Decay)", 
                 "Edwards-indvD-HS", 
                 "Edwards-D-HS", 
                 "Edwards-D-HS\n(No Synergy)", 
                 "Gompertz-D-HS",
                 "Gompertz-D-HS\n(No Synergy)")

model_colours <- c("#b2df8a",
                   "#fb9a99",
                   "#1f78b4",
                   "#6a3d9a",
                   "#a6cee3",
                   "#cab2d6", 
                   "#3690c0",
                   "#8c6bb1", 
                   "#f768a1",
                   "#dd3497",
                   "#ae017e",
                   "#88419d",
                   "#810f7c",
                   "#4d004b",
                   "#6a51a3",
                   "#045a8d",
                   "#045a8d")

#
#  Read in model CV performance
#
drug_free_file_path <- "output/drug_free/"
indv_drug_file_path <- "output/indv_drug/"
all_drug_file_path <- "output/all_drug/"
file_paths <- list(drug_free_file_path, indv_drug_file_path, all_drug_file_path)
fit_metrics <- list()
for (i in 1:length(file_paths)) {
  list_fit_metrics <- lapply(model_lists[[i]], function(x){readRDS(file=paste(file_paths[[i]], x, "_fit_stats.rds", sep = "")) %>% add_column(model = x) %>% add_column(condition = conditions[i])})
  fit_metrics[[i]] <- do.call(rbind, list_fit_metrics)
}
total_fit_metrics <- do.call(rbind, fit_metrics)
#
# Plot predictive performance
#
total_fit_metrics %>%
  filter(testing==TRUE) %>%
  group_by(metric, model, condition) %>% 
  mutate(mean_mu = mean(mean)) %>%   # calc. mean and se over folds
  mutate(mean_sem = sd(mean)/sqrt(3)) %>%
  ungroup() %>%
  unique() %>%
  filter(!model%in%c("gompertz_drug_hs_reg", "gompertz_drug_no_syn_hs_reg")) %>%  # filter out models only shown in supp.
  filter(!grepl("decay", model)) %>%
  mutate(model = factor(model, levels=unlist(model_lists), labels = model_names)) %>%
  mutate(condition = factor(condition, levels=conditions, labels=c("Drug Free", "Excluding Combinations", "All"))) %>%
  filter(metric%in%c("RMSE", "LPD")) %>%
  group_by(metric) %>%
  mutate(diff = max(mean)-mean_mu + 1) %>%
  mutate(diff_fold = max(mean)-mean + 1) %>%
  ungroup() %>%
  mutate(values = case_when(metric=="LPD" ~ diff, metric=="RMSE" ~ mean_mu)) %>%
  mutate(values_fold = case_when(metric=="LPD" ~ diff_fold, metric=="RMSE" ~ mean)) %>%
  mutate(metric = case_when(metric == "LPD" ~ "Relative LPD", .default=metric)) %>%
  mutate(metric = factor(metric, levels = c("Relative LPD", "RMSE"), labels = c(expression("log"[10]~"Relative LPD"), expression("log"[10]~"RMSE")))) %>%
  ggplot(aes(x = model, y = values, color=model)) +
    geom_errorbar(aes(ymin=values-mean_sem, ymax=values+mean_sem), width=0.1, colour="#737373", alpha=0.9, size=0.5) +
    geom_point(size=2) +
    geom_jitter(aes(y = values_fold, shape=factor(fold)), size=1.5, width=0.25, alpha=0.7) +
    ggh4x::facet_grid2(metric ~ condition, scales = "free", independent = "y", switch="y", labeller = labeller(metric = label_parsed)) +
    theme_bw(base_size=11) +
    theme(strip.placement = "outside", strip.background = element_blank(), axis.text.x = element_text(angle=90, vjust=0.7, hjust=0.7), strip.text.x=element_text(size=11), strip.text.y=element_text(size=11), legend.position="bottom") +
    guides(color = "none") +
    xlab("Models") +
    ylab("") +
    scale_color_manual(breaks = model_names, values = model_colours) +
    scale_shape_manual(name="Fold", breaks = 1:3, values = c(18, 17, 15)) +
    scale_y_continuous(breaks=c(10^(-1.8), 10^(-1.75), 10^(-1.6), 10^(-1.5), 1e-1, 10^(-0.5), 1, 1.25, 1.5, 2, 5), trans="log10",  labels = c(-1.8, -1.75, -1.6, -1.5, -1, -0.5, 0, 0.097, 0.176, 0.301, 0.699)) +
    force_panelsizes(cols = c(3, 3, 1)) -> p
tiff("figures/fig2.tif", width = 20, height = 17, units = "cm", res=300)
p
dev.off()
#
# Plot fitting and predictive performance of all models
#
total_fit_metrics %>%
  group_by(metric, model, condition, testing) %>%
  mutate(mean_mu = mean(mean)) %>%   # calc. mean and se over folds
  mutate(mean_sem = sd(mean)/sqrt(3)) %>%
  ungroup() %>%
  unique() %>%
  mutate(model = factor(model, levels=unlist(model_lists), labels = model_names)) %>%
  mutate(condition = factor(condition, levels=conditions, labels=c("Drug Free", "Excluding Combinations", "All"))) %>%
  filter(metric%in%c("RMSE", "LPD")) %>%
  group_by(metric, testing) %>%
  mutate(diff = max(mean)-mean_mu + 1) %>%
  mutate(diff_fold = max(mean)-mean + 1) %>%
  ungroup() %>%
  mutate(values = case_when(metric=="LPD" ~ diff, metric=="RMSE" ~ mean_mu)) %>%
  mutate(values_fold = case_when(metric=="LPD" ~ diff_fold, metric=="RMSE" ~ mean)) %>%
  mutate(metric = case_when(metric == "LPD" ~ "Relative LPD", .default=metric)) %>%
  mutate(metric_testing = case_when(metric == "RMSE" & !testing ~ "Train RMSE",
                                    metric == "RMSE" & testing ~ "Test RMSE",
                                    metric == "Relative LPD" & testing ~ "Test Relative LPD")) %>%
  mutate(metric_testing = factor(metric_testing, levels = c("Test Relative LPD", "Test RMSE", "Train RMSE"), labels = c(expression("Test log"[10]~"Relative LPD"), expression("Test log"[10]~"RMSE"), expression("Train log"[10]~"RMSE")))) %>%
  ggplot(aes(x = model, y = values, color=model)) +
    geom_errorbar(aes(ymin=values-mean_sem, ymax=values+mean_sem), width=0.1, colour="#737373", alpha=0.9, size=0.5) +
    geom_point(size=2) +
    geom_jitter(aes(y = values_fold, shape=factor(fold)), size=1.5, width=0.25, alpha=0.7) +
    ggh4x::facet_grid2(metric_testing ~ condition, scales = "free", independent = "y", switch="y", labeller = labeller(metric_testing = label_parsed)) +
    theme_bw(base_size=11) +
    theme(strip.placement = "outside", strip.background = element_blank(), axis.text.x = element_text(angle=90, vjust=0.7, hjust=0.7), strip.text.x=element_text(size=11), strip.text.y=element_text(size=11), legend.position="bottom") +
    guides(color = "none") +
    xlab("Models") +
    ylab("") +
    scale_color_manual(breaks = model_names, values = model_colours) +
    scale_shape_manual(name="Fold", breaks = 1:3, values = c(18, 17, 15)) +
    scale_y_continuous(breaks=c(10^(-2.05), 1e-2, 10^(-1.95), 10^(-1.9), 10^(-1.8), 10^(-1.75), 10^(-1.6), 10^(-1.5), 1e-1, 10^(-0.5), 1, 1.25, 1.5, 2, 5), trans="log10",  labels = c(-2.05, -2, -1.95, -1.9, -1.8, -1.75, -1.6, -1.5, -1, -0.5, 0, 0.097, 0.176, 0.301, 0.699)) +
    force_panelsizes(cols = c(6, 7, 4)) -> p2

tiff("figures/figs1.tif", width = 20, height = 22, units = "cm", res=300)
p2
dev.off()
