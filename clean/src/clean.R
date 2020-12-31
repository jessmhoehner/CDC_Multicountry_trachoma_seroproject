#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JH
# Maintainer(s): JH
# License: 2020, GPL v2 or later
#
# -----------------------------------------------------------
# multicountrytrachseroproject/clean/src/clean.R

## Cleaning code for                                               ##
## "Comparison of multiple tests for determination of              ##
## seroconversion rates to the Chlamydia trachomatis antigen       ##
## Pgp3: a multi-country analysis"                                 ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##   
## but please maintain appropriate accreditation)                  ##
##                                                                 ##   
## Jessica Hoehner (nee Randall)                                             ##
## github.com/jessmhoehner                                     ##

## this script reads in individual .csvs with age, titre, and sero_pos columns
## after reading in it make sure none of the sheets are empty, they each have
## 3 columns 

## next, the cleaned dataframes in each sheet are used to calculate prevalence
## proportions for each age group and 95% binomial confidence estimates
## sheets containing the prevalences and CI estimates are exported to the plot
## task's input folder for plotting age seroprevalence curves


# load libraries 
pacman::p_load("here","readr", "janitor", "tidyverse", 
               "assertr", "purrr", "binom")

# specify file structure of inputs and outputs

files <- list(
  drc1_Ct694 = here("clean/input/DRC1Ct694.csv"), 
  drc1_LFA = here("clean/input/DRC1LFA.csv"),
  drc1_MBA = here("clean/input/DRC1MBA.csv"),
  drc2_Ct694 = here("clean/input/DRC2Ct694.csv"),
  drc2_LFA = here("clean/input/DRC2LFA.csv"),
  drc2_MBA = here("clean/input/DRC2MBA.csv"),
  togoLFAf_41 = here("clean/input/TogoLFAfield_40001.csv"), 
  togoLFAf_42 = here("clean/input/TogoLFAfield_40002.csv"),
  togoLFAg_41 = here("clean/input/TogoLFAgold_40001.csv"), 
  togoLFAg_42 = here("clean/input/TogoLFAgold_40002.csv"), 
  togoLFAl_41 = here("clean/input/TogoLFAlatex_40001.csv"), 
  togoLFAl_42 = here("clean/input/TogoLFAlatex_40002.csv"), 
  togoMBAc_41 = here("clean/input/TogoMBAct694_40001.csv"), 
  togoMBAc_42 = here("clean/input/TogoMBAct694_40002.csv"), 
  togoMBAp_41 =  here("clean/input/TogoMBAPgp3_40001.csv"), 
  togoMBAp_42 = here("clean/input/TogoMBAPgp3_40002.csv")

)

stopifnot(is_empty(files) != TRUE & length(files) == 16)
## Read in data

## creates a list of all files as connections
fileslist <- list(files$drc1_Ct694 , files$drc1_LFA, files$drc1_MBA, 
                  files$drc2_Ct694, files$drc2_LFA, files$drc2_MBA, 
                  files$togoLFAf_41, files$togoLFAf_42, files$togoLFAg_41, 
                  files$togoLFAg_42, files$togoLFAl_41, files$togoLFAl_42, 
                  files$togoMBAc_41, files$togoMBAc_42, files$togoMBAp_41, 
                  files$togoMBAp_42)

stopifnot(length(fileslist) == 16)

# iterates over list of files, cleans the names of the columns, checks for 
# only 3 columns in each file, and that no values are missing
cleanlist <- lapply(fileslist, function(x) {
  
  x_df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
  clean_names()
  
  x_df  %>%
    verify(ncol(x_df) == 3) %>%
    verify(is.na(x_df) == FALSE)
  
})

stopifnot(length(cleanlist) == 16)

# add unique names to each df for easy export later on

df_names <- c("drc1_Ct694", "drc1_LFA", "drc1_MBA", "drc2_Ct694", "drc2_LFA", 
              "drc2_MBA", "togoLFAf_41", "togoLFAf_42", "togoLFAg_41", 
              "togoLFAg_42", "togoLFAl_41", "togoLFAl_42", "togoMBAc_41", 
              "togoMBAc_42", "togoMBAp_41", "togoMBAp_42")

names(cleanlist) <- df_names

## using cleanlist, we extract each df and save and export result 
## to the model and the observed tasks respectively##

#start i loop
for (i in seq_along(cleanlist)) {
  
  df <- as.data.frame(pluck(cleanlist, i))
  
  write_excel_csv(df, quote = FALSE, 
                  path = 
                      here(paste("model/input/",names(cleanlist)[i],
                                 "_cleanmod.csv", sep = "")))
  
  write_excel_csv(df, quote = FALSE, 
                    path = 
                      here(paste("observed/input/",names(cleanlist)[i],
                                 "_cleanobs.csv", sep = "")))
  
  write_excel_csv(df, quote = FALSE, 
                  path = 
                    here(paste("plot/input/",names(cleanlist)[i],
                               "_cleanobs.csv", sep = "")))
  
  #message to let the user know that each iteration has completed
  print(paste0("Cleaning for dataset ",names(cleanlist)[i],
               " has completed successfully."))
  
  } # close i loop

# done
