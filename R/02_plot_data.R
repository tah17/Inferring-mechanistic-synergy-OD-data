rm(list = ls())
library(tidyverse)
seed <- 404806
set.seed(seed)
#
# reads in processed data
#
readRDS(file = "data/OD.Rda") %>%
  drop_na(drug_conds) %>%  # drop wells with nothing in them
  mutate(drug_conds = factor(drug_conds, levels = c("2X RPMI + dH2O", 
                                                    "RPMI", 
                                                    "8 AFG", 
                                                    "0.03 MGX", 
                                                    "0.25 5FC", 
                                                    "8 AFG +\r\n0.03 MGX", 
                                                    "8 AFG +\r\n0.25 5FC"))) %>%
  select(-c(strain)) %>%
  separate(drug_name, c("drug1_name", "drug2_name"), "_", remove=TRUE) %>%
  pivot_longer(c(drug1_name, drug2_name), names_to="drug", values_to = "drug_name") %>% 
  mutate(drug_name = factor(drug_name, levels = c(NA, "AFG", "MGX", "5FC"), labels = c("None", "AFG", "MGX", "5FC"), exclude=NULL)) %>%
  pivot_wider(names_from = "drug", values_from = "drug_name") -> data
#
# plot data
#
pal_fill <- scales::brewer_pal(palette = "Dark2")(5)  # colour palette for drug conds
data %>%
  unite("drug_name", c("drug1_name", "drug2_name"), remove=TRUE) %>%
  mutate(drug_name = case_when(drug_conds=="2X RPMI + dH2O" ~ "Blanks",
                               drug_conds=="RPMI" ~ "RPMI-only", 
                               drug_name=="AFG_5FC" ~ "AFG+5FC",
                               drug_name=="AFG_MGX" ~ "AFG+MGX",
                               grepl("_None", drug_name) ~ str_replace(drug_name, "_None", ""),
                               .default=drug_name)) %>%
  mutate(drug_name = factor(drug_name, levels = c("Blanks", "RPMI-only", "AFG", "MGX", "5FC", "AFG+MGX", "AFG+5FC"))) %>%
  arrange(drug_name) %>%
  ggplot(aes(x = time, y = OD,  group = well, color = drug_name)) +
  facet_wrap(. ~ drug_name, nrow=2) +
  geom_point(size=0.5) +
  geom_line(alpha=0.5) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.placement = "outside") +
  ylab("OD") +
  xlab("Time [hrs]") +
  scale_colour_manual(name="Drug Condition", values=c("#bdbdbd", "#737373", pal_fill)) -> p

tiff("figures/fig1.tif", width = 22, height = 12, units = "cm", res=300)
p
dev.off()

