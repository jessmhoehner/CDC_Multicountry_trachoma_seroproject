#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License: 2020, EICC, GPL v2 or later
#
# -----------------------------------------------------------
# multicountrytrachseroproject/adjust/src/adjust.R

## Age adjustment code for                                         ##
## "Comparison of multiple tests for determination of              ##
## seroconversioni rates to the Chlamydia trachomatis antigen      ##
## Pgp3: a multi-country analysis"                                 ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##   
## but please maintain appropriate accreditation)                  ##
##                                                                 ##   
## Jessica Randall                                                 ##
## github.com/jessmrandall                                         ##

# load libraries 
pacman::p_load("here","readr", "janitor", "tidyverse", "assertr")

# specify file structure of inputs and outputs

files <- list(drc1_Ct694_prevs = 
                here("/adjust/input/DRC1Ct694_obs_adjust_df.csv"),
              ageweights_drc_clean = 
                here("/adjust/output/age_weights_DRC.csv"),
              age_adj_prevs = 
                here("/adjust/output/DRCCt694_ageadjustedprevs"),
              values_table = 
                here("/adjust/output/DRCCt694_overallvaluestable")
)

# read in observed prevalence data

drc1_ct694_df <- read_csv(files$drc1_Ct694_prevs, col_names = TRUE, na = "NA") %>%
  select(-c("age_bins_mid"))

# load age weights obtained from Tropical Data###################

age	<- as.numeric(c(1,2,3,4,5,6,7,8,9))

weights <- as.numeric(c(0.11823272, 0.11823272, 0.11823272, 0.11823272, 
                        0.10541382, 0.10541382, 0.10541382, 0.10541382, 
                        0.10541382))

age_weights_drc <- data.frame(age, weights) %>%
  write_excel_csv(files$ageweights_drc_clean)

##########################################################
# adjust data by age weights and export

age_adj_obs_data <- drc1_ct694_df %>%
  mutate(aa_prev = ((med*100)*(age_weights_drc$weights)), 
         aa_low95 = ((low_95*100)*(age_weights_drc$weights)), 
         aa_high95 = ((high_95*100)*(age_weights_drc$weights)))

write_excel_csv(age_adj_obs_data, files$drc1_ct694_aa_obsdf)

##########################################################
#calculate and export table of overall raw and adjusted values

# raw #
# overall raw % pos
raw_ppos <- round(mean((age_adj_obs_data$med)*100), digits = 2)
print(paste0("The overall mean raw percentage of positive tests by Ct694 in Manono, DRC (40586) is ", raw_ppos,"%" ))

# overall raw lower CI mean #
raw_plcl <- round(mean((age_adj_obs_data$low_95*100)), digits = 2)
print(paste0("The overall mean raw lower 95% confidence interval estimate around the percentage of positive tests by Ct694 in Manono, DRC (40586) is ", raw_plcl,"%" ))

# overall raw upper CI mean #
raw_pucl <- round(mean((age_adj_obs_data$high_95*100)), digits = 2)
print(paste0("The overall mean raw upper 95% confidence interval estimate around the percentage of positive tests by Ct694 in Manono, DRC (40586) is ", raw_pucl,"%" ))

#age adjusted#
# overall age adjusted %pos 
aa_ppos <- round(sum(age_adj_obs_data$aa_prev), digits = 2)
print(paste0("The overall sum age-adjusted percentage of positive tests by Ct694 in Manono, DRC (40586) is ", raw_ppos,"%" ))

# overall age adjusted lower CI sum #
aa_plcl <- round(sum(age_adj_obs_data$aa_low95), digits = 2)
print(paste0("The overall sum age-adjusted lower 95% confidence interval estimate around the percentage of positive tests by Ct694 in Manono, DRC (40586) is ", aa_plcl,"%" ))

# overall age adjusted upper CI sum #
aa_pucl <- round(sum(age_adj_obs_data$aa_high95), digits = 2)
print(paste0("The overall sum age-adjusted lower 95% confidence interval estimate around the percentage of positive tests by Ct694 in Manono, DRC (40586) is ", aa_pucl,"%" ))

#create a table of all summary values
overall <- data.frame(raw_ppos, raw_plcl, raw_pucl, 
                     aa_ppos, aa_plcl, aa_pucl) %>%
  transmute(percent_pos_raw = as.numeric(raw_ppos), 
         percent_lcl_raw = as.numeric(raw_plcl), 
         percent_ucl_raw = as.numeric(raw_pucl),
         percent_pos_adj = as.numeric(aa_ppos), 
         percent_lcl_adj = as.numeric(aa_plcl), 
         percent_ucl_adj = as.numeric(aa_pucl))

# export table of raw and adjusted overall values
write_excel_csv(overall, files$values_table)

#done