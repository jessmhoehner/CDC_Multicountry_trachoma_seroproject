!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License: 2020, EICC, GPL v2 or later
#
# -----------------------------------------------------------
# NCV19/clean/src/clean.R

## Cleaning code for                                               ##
## "Comparison of multiple tests for determination of              ##
## seroconversion rates to the Chlamydia trachomatis antigen       ##
## Pgp3: a multi-country analysis"                                 ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##   
## but please maintain appropriate accreditation)                  ##
##                                                                 ##   
## Jessica Randall                                                 ##
## github.com/jessmrandall                                         ##


# load libraries 
pacman::p_load("here","readr", "janitor", 
               "tidyverse", "assertr", "purrr", 
               "binom")

# specify file structure of inputs and outputs

files <- list(
  drc1_Ct694 = here("/clean/input/DRC1Ct694.csv"), 
  drc1_LFA = here("/clean/input/DRC1LFA.csv"),
  drc1_MBA = here("/clean/input/DRC1MBA.csv"),
  drc2_Ct694 = here("/clean/input/DRC1Ct694.csv"), 
  drc2_LFA = here("/clean/input/DRC1LFA.csv"),
  drc2_MBA = here("/clean/input/DRC1MBA.csv"),
  
  togoLFAf_41 = here("/clean/input/TogoLFAfield_40001.csv"), 
  togoLFAf_42 = here("/clean/input/TogoLFAfield_40002.csv"), 
  togoLFAg_41 = here("/clean/input/TogoLFAgold_40001.csv"), 
  togoLFAg_42 = here("/clean/input/TogoLFAgold_40002.csv"), 
  togoLFAl_41 = here("/clean/input/TogoLFAlatex_40001.csv"), 
  togoLFAl_42 = here("/clean/input/TogoLFAlatex_40002.csv"), 
  togoMBAc_41 = here("/clean/input/TogoMBAct694_40001.csv"), 
  togoMBAc_42 = here("/clean/input/TogoMBAct694_40002.csv"), 
  togoMBAp_41 =  here("/clean/input/TogoMBAct694_40001.csv"), 
  togoMBAp_42 = here("/clean/input/TogoMBAct694_40002.csv"),
  
  
  drc1_Ct694_clean = here("/model/input/DRC1Ct694_clean.csv"),
  drc1_ct694_observed = here("/plot/input/DRC1Ct694_observed_df.csv"),
  
  drc1_LFA_clean = here("/model/input/DRC1LFA_clean.csv"),
  
  
  
  drc1_MBA_clean = here("/model/input/DRC1MBA_clean.csv"),
  
  
  
  drc2_Ct694_clean = here("/model/input/DRC1Ct694_clean.csv"),
  
  
  
  drc2_LFA_clean = here("/model/input/DRC1LFA_clean.csv"),
  
  
  
  drc2_MBA_clean = here("/model/input/DRC1MBA_clean.csv"),
  
  
  
  
  togoLFAf_41_clean = here("/model/input/TogoLFAfield_40001_clean.csv"),
  
  
  
  togoLFAf_42_clean = here("/model/input/TogoLFAfield_40002_clean.csv"),
  
  
  
  togoLFAg_41_clean = here("/model/input/TogoLFAgold_40001_clean.csv"), 
  
  
  
  togoLFAg_42_clean = here("/model/input/TogoLFAgold_40002_clean.csv"), 
  
  
  
  togoLFAl_41_clean = here("/model/input/TogoLFAlatex_40001_clean.csv"), 
  
  
  
  togoLFAl_42_clean = here("/model/input/TogoLFAlatex_40002_clean.csv"),
  
  
  
  togoMBAc_41_clean = here("/model/input/TogoMBAct694_40001_clean.csv"), 
  
  
  
  togoMBAc_42_clean = here("/model/input/TogoMBAct694_40002_clean.csv"),
  
  
  
  togoMBAp_41_clean = here("/model/input/TogoMBAct694_40001_clean.csv"),
  
  
  
  togoMBAp_42_clean = here("/model/input/TogoMBAct694_40002_clean.csv")
)

stopifnot(is_empty(files) != TRUE & length(files) == 33)
## Read in data

# DRC

# DRC1_CT694

drc1_ct694_df <- read_csv(files$drc1_Ct694, col_names = TRUE, na = "NA") %>%
  clean_names()

drc1_ct694_df %>%
  verify(ncol(drc1_ct694_df) == 3 & nrow(drc1_ct694_df) == 1496) %>%
  verify(is.na(drc1_ct694_df) == FALSE) %>%
  write_excel_csv(files$drc1_Ct694_clean)

# add in age weighting here, look in "fitted_model_output for DRC1" on the 
# bottom left hand side (~ line 18)

# age bin the data for plotting
age_bins <- seq(from=0, to=10, by=1)
age_bins_mid <- seq(from=0.5, to=9.5, by=1) 

N_bins <- length(age_bins) - 1 

# initialize empty dataframe to fill with binomial confidence intervals
sp_bins <- data.frame(med = numeric(0), 
                      low_95 = numeric(0), 
                      high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix					  
for(i in 1:N_bins)
{
  index <- which(drc1_ct694_df$age > age_bins[i] & 
                   drc1_ct694_df$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc1_ct694_df[index,3]
  
  sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                    nrow(temp_AB), 
                                                    method="wilson", 
                                                    seed = seed)[1,4:6]))
}

stopifnot(is_empty(sp_bins) == FALSE)
stopifnot(nrow(sp_bins) == 9 & ncol(sp_bins) == 3)

#setup data frames for plots and export
sp_bins_df <- sp_bins %>%
  mutate(age = as.factor(row.names(sp_bins)))

age_bins_mid_df <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(sp_bins_df$age))

observed_data_df <- left_join(sp_bins_df, age_bins_mid_df, by = "age") %>%
  mutate(age = as.numeric(age))

write_excel_csv(observed_data_df, files$drc1_ct694_observed)

# drc1_LFA 
# drc1_MBA 
# drc2_Ct694 
# drc2_LFA 
# drc2_MBA 

# Togo

# done


# drc1_LFA 
# drc1_MBA 
# drc2_Ct694 
# drc2_LFA 
# drc2_MBA 

# Togo

# done