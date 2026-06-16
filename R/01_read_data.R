#
# Script that reads in the spreadsheet data from the OD reader.
#
#
rm(list = ls())
library(tidyverse)
library(readxl)
library(reshape2)
library(janitor)
library(RColorBrewer)
library(ggrepel)

# Read in data set 1 -------------------------------------------------
excel_file <- "data/raw_OD.xlsx"
sheet_name <- "Plate\ Layout"
data_range <- "C9:N17"
strain <- tibble(strain = rep("B12663", 18), well = unlist(lapply(3:5, function(x) paste(LETTERS[x], 4:9, sep = ""))))
#
# creates meta-data
#
read_excel(excel_file, sheet = sheet_name, range=data_range) %>%
  mutate(well = LETTERS[row_number()]) %>%
  pivot_longer(!well) %>%
  unite("well", c(well, name), sep = "") %>%
  rename(drug_conds = value) %>%
  mutate(drug_conds = case_when(well%in%c(paste("A", 4:9, sep=""), unlist(lapply(c(2, 6:7), function(x) paste(LETTERS[x], 2:11, sep = ""))), unlist(lapply(3:5, function(x) paste(LETTERS[x], c(2:3, 10:11), sep = "")))) ~ "2X RPMI + dH2O", .default = drug_conds)) %>%
  mutate(blanks = case_when(drug_conds%in%c("2X RPMI + dH2O", NA) ~ TRUE, !(drug_conds%in%c("2X RPMI + dH2O", NA)) ~ FALSE)) %>%
  full_join(strain) -> meta_data
#
# reads in OD
#
data_range <- "A33:AX130"
sheet_name <- "Result\ sheet"
df <- read_excel(excel_file, sheet = sheet_name, range = data_range) %>% 
  filter(row_number() != 1) %>%
  rename(well = `Time [s]`) %>%
  pivot_longer(!well, names_to = "time_s", values_to = "OD") %>%
  mutate(time_s = as.numeric(time_s))  %>%
  mutate(time = time_s/(60*60)) %>%  # converts time to hours
  select(-time_s)
#
#  combines meta-data with OD
#
raw_data <- full_join(meta_data, df)
#
#  quickly plots data
#
raw_data %>%
  ggplot(aes(x = time, y = OD, group = well, color = drug_conds)) +
  geom_point(size=0.5) +
  geom_line(alpha=.5) +
  facet_wrap(blanks ~ .) +
  # scale_y_log10() +
  theme_bw(base_size = 20)

# Process data sets -------------------------------------------------
#
# processes data to get drug concs
#
raw_data %>%
  mutate(drug_conds = case_when(drug_conds=="GC" ~ "RPMI", .default = drug_conds)) %>%
  separate(drug_conds, c("drug1_conc", "drug1_name_drug2"), " ", remove = FALSE, extra="merge") %>%
  separate(drug1_name_drug2, c("drug1_name", "drug2"), " \\+ | \\+\r\n", extra="merge") %>% 
  separate(drug2, c("drug2_conc", "drug2_name"), " ", extra="merge") %>%
  mutate(drug1_conc = case_when(!(drug_conds%in%c("2X RPMI + dH2O", "RPMI", NA)) ~ as.numeric(drug1_conc), .default=NA)) %>%
  mutate(drug2_conc = case_when(!(drug_conds%in%c("2X RPMI + dH2O", "RPMI", NA)) ~ as.numeric(drug2_conc), .default=NA)) %>%
  mutate(drug1_name = case_when(drug_conds%in%c("2X RPMI + dH2O", "RPMI", NA) ~ NA, .default=drug1_name)) %>%
  mutate(drug2_name = case_when(drug_conds%in%c("2X RPMI + dH2O", "RPMI", NA) ~ NA, .default=drug2_name)) %>%
  unite(drug_name, c(drug1_name, drug2_name), na.rm=TRUE) %>%
  unite(drug_conc, c(drug1_conc, drug2_conc), na.rm=TRUE) %>%
  mutate(drug_name = case_when(drug_name=="" ~ NA, .default=drug_name)) %>%
  mutate(drug_conc = case_when(drug_conc=="" ~ NA, .default=drug_conc)) -> data
#
# saves processed data
#
saveRDS(data, file="data/OD.Rda")

