#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JH
# Maintainer(s): JH
# License: 2020, GPL v2 or later
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
## Jessica Hoehner                                                 ##
## github.com/jessmhoehner                                         ##

## this script adjusts the observed sero prevalences in the input  ##
## csvs by age weights provided from Tropical Data                 ##

# load libraries 
pacman::p_load("here","readr", "janitor", "tidyverse", "assertr")

# specify file structure of inputs and outputs
# have Togo data in the input folder but no weights for Togo specifically

files <- list(drc1_Ct694_prevs = 
                here("adjust/input/drc1_Ct694_obs_df.csv"),
              drc1_Ct694_prevs = here("adjust/input/DRC1Ct694.csv"), 
              drc1_LFA_prevs = here("adjust/input/drc1_LFA_obs_df.csv"),
              drc1_MBA_prevs = here("adjust/input/drc1_MBA_obs_df.csv"),
              drc2_Ct694_prevs = here("adjust/input/drc2_Ct694_obs_df.csv"),
              drc2_LFA_prevs = here("adjust/input/drc2_LFA_obs_df.csv"),
              drc2_MBA_prevs = here("adjust/input/drc2_MBA_obs_df.csv"),
              
              ageweights_drc_clean = 
                here("adjust/output/age_weights_DRC.csv"),
              age_adj_prevs = 
                here("adjust/output/DRCCt694_ageadjustedprevs.csv"),
              values_table = 
                here("adjust/output/DRCCt694_overallvaluestable.csv")
)

# read in observed prevalence data
# currently only have weights for DRC

## creates a list of all files as connections
fileslist <- list(files$drc1_Ct694_prevs, files$drc1_LFA_prevs, 
                  files$drc1_MBA_prevs,files$drc2_Ct694_prevs, 
                  files$drc2_LFA_prevs,files$drc2_MBA_prevs)

stopifnot(length(fileslist) == 6)

# iterates over list of files, cleans the names of the columns, checks for 
# only 5 columns in each file, and that no values are missing
adjlist <- lapply(fileslist, function(x) {
  
  x_df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
    clean_names()
  
  x_df  %>%
    verify(ncol(x_df) == 5) %>%
    verify(is.na(x_df) == FALSE)
  
})

stopifnot(length(adjlist ) == 6)

# add unique names to each df for easy export later on

df_names <- c("drc1_CT694", 
              "drc1_LFA", 
              "drc1_PgP3", 
              "drc2_CT694", 
              "drc2_LFA", 
              "drc2_Pgp3")

names(adjlist) <- df_names

## using cleanlist, we extract each df and save and export result 
## to the model and the observed tasks respectively##

#start i loop
for (i in seq_along(adjlist)) {
  
  #messages for the user to keep them aware of model progress
  start_time <- Sys.time()
  print(paste0("Age weight adjustment for dataset ",names(adjlist)[i]," has now begun..."))
  
  df <- as.data.frame(pluck(adjlist, i))
  
  # load age weights obtained from Tropical Data###################
  
  age	<- as.numeric(c(1,2,3,4,5,6,7,8,9))
  
  #manually added weights as list
  weights_drc <- as.numeric(c(0.11823272, 0.11823272, 0.11823272, 0.11823272, 
                          0.10541382, 0.10541382, 0.10541382, 0.10541382, 
                          0.10541382))
  
  #xport weights to output folder
  age_weights_drc <- data.frame(age, weights_drc)
  
  write_excel_csv(age_weights_drc, quote = FALSE, 
                  path = 
                      here(paste("adjust/output/",names(adjlist)[i],"_ageweights_df.csv", sep = "")))

  ##########################################################
  
  #export weights as a dataframe for each test type and community
  age_adj_obs_data <- df %>%
    mutate(aa_prev = ((med*100)*(age_weights_drc$weights)), 
           aa_low95 = ((low_95*100)*(age_weights_drc$weights)), 
           aa_high95 = ((high_95*100)*(age_weights_drc$weights)))

  write_excel_csv(age_adj_obs_data, quote = FALSE, 
                  path = 
                    here(paste("adjust/output/",names(adjlist)[i],"_adjusted_df.csv", sep = "")))

  #calculate and export table of overall raw and adjusted values
  
  # raw #
  # overall raw % pos
  raw_ppos <- round(mean((age_adj_obs_data$med)*100), digits = 2)
  
  # overall raw lower CI mean #
  raw_plcl <- round(mean((age_adj_obs_data$low_95*100)), digits = 2)
  
  # overall raw upper CI mean #
  raw_pucl <- round(mean((age_adj_obs_data$high_95*100)), digits = 2)
  
  #age adjusted#
  # overall age adjusted %pos 
  aa_ppos <- round(sum(age_adj_obs_data$aa_prev), digits = 2)
  
  # overall age adjusted lower CI sum #
  aa_plcl <- round(sum(age_adj_obs_data$aa_low95), digits = 2)
  
  # overall age adjusted upper CI sum #
  aa_pucl <- round(sum(age_adj_obs_data$aa_high95), digits = 2)
  
  #create a table of all summary values
  summarytable <- data.frame(raw_ppos, raw_plcl, raw_pucl, 
                       aa_ppos, aa_plcl, aa_pucl) %>%
    transmute(percent_pos_raw = as.numeric(raw_ppos), 
           percent_lcl_raw = as.numeric(raw_plcl), 
           percent_ucl_raw = as.numeric(raw_pucl),
           percent_pos_adj = as.numeric(aa_ppos), 
           percent_lcl_adj = as.numeric(aa_plcl), 
           percent_ucl_adj = as.numeric(aa_pucl))
  
  # export table of raw and adjusted overall values
  write_excel_csv(summarytable, quote = FALSE, 
                  path = 
                    here(paste("adjust/output/",names(adjlist)[i],"_summarytable.csv", sep = "")))

#message to let the user know that each iteration has completed
print(paste0("Adjustment for dataset ",names(adjlist)[i]," has completed successfully."))

}

#done
