#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License: 2020, EICC, GPL v2 or later
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
## Jessica Randall                                                 ##
## github.com/jessmrandall                                         ##

## this script reads in individual .csvs with age, titre, and sero_pos columns
## after reading in it make sure none of the sheets are empty, they each have
## 3 columns and however many columns they are expected to have (tested with 
## unit tests), and they are exported to the model task input folder

## next, the cleaned dataframes in each sheet are used to calculate prevalence
## proportions for each age group and 95% binomial confidence estimates
## sheets containing the prevalences and CI estimates are exported to the plot
## task's input folder for plotting age seroprevalence curves


# load libraries 
pacman::p_load("here","readr", "janitor", "tidyverse", 
               "assertr", "purrr", "binom")

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
  
  drc1_Ct694_clean = here("/observed/input/DRC1Ct694_clean.csv"),
  drc1_LFA_clean = here("/observed/input/DRC1LFA_clean.csv"),
  drc1_MBA_clean = here("/observed/input/DRC1MBA_clean.csv"),
  drc2_Ct694_clean = here("/observed/input/DRC2Ct694_clean.csv"),
  drc2_LFA_clean = here("/observed/input/DRC2LFA_clean.csv"),
  drc2_MBA_clean = here("/observed/input/DRC2MBA_clean.csv"),
  togoLFAf_41_clean = here("/observed/input/TogoLFAfield_40001_clean.csv"),
  togoLFAf_42_clean = here("/observed/input/TogoLFAfield_40002_clean.csv"),
  togoLFAg_41_clean = here("/observed/input/TogoLFAgold_40001_clean.csv"),
  togoLFAg_42_clean = here("/observed/input/TogoLFAgold_40002_clean.csv"), 
  togoLFAl_41_clean = here("/observed/input/TogoLFAlatex_40001_clean.csv"), 
  togoLFAl_42_clean = here("/observed/input/TogoLFAlatex_40002_clean.csv"),
  togoMBAc_41_clean = here("/observed/input/TogoMBAct694_40001_clean.csv"),
  togoMBAc_42_clean = here("/observed/input/TogoMBAct694_40002_clean.csv"),
  togoMBAp_41_clean = here("/observed/input/TogoMBAct694_40001_clean.csv"),
  togoMBAp_42_clean = here("/observed/input/TogoMBAct694_40002_clean.csv")

)

stopifnot(is_empty(files) != TRUE & length(files) == 32)
## Read in data

#set seed for reproducibility
seed <- 22315

## creates a list of all files
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

#check number of rows (vary by data frame) and export to next task as 
# individually named dataframe    #

### DRC ###
###
drc1_ct694 <- as.data.frame(cleanlist[[1]]) 

drc1_ct694  %>%
  verify(nrow(drc1_ct694) == 1496) %>%
  write_excel_csv(files$drc1_Ct694_clean)

### 
drc1_LFA <- as.data.frame(cleanlist[[2]]) 

drc1_LFA %>%
  verify(nrow(drc1_LFA) == 1494) %>%
  write_excel_csv(files$drc1_LFA_clean)

###

drc1_MBA <- as.data.frame(cleanlist[[3]]) 

drc1_MBA %>%
  verify(nrow(drc1_MBA) == 1496) %>%
  write_excel_csv(files$drc1_MBA_clean)

###
drc2_Ct694 <- as.data.frame(cleanlist[[4]])

drc2_Ct694 %>%
  verify(nrow(drc2_Ct694) == 1496) %>%
  write_excel_csv(files$drc2_Ct694_clean)

###

drc2_LFA <- as.data.frame(cleanlist[[5]]) 

drc2_LFA %>%
  verify(nrow(drc2_LFA) == 1494) %>%
  write_excel_csv(files$drc2_LFA_clean)

###

drc2_MBA <- as.data.frame(cleanlist[[6]]) 

drc2_MBA %>%
  verify(nrow(drc2_MBA) == 1496) %>%
  write_excel_csv(files$drc2_MBA_clean)
###

###Togo###
###

togoLFAf41 <- as.data.frame(cleanlist[[7]]) 

togoLFAf41 %>%
  verify(nrow(togoLFAf41) == 972) %>%
  write_excel_csv(files$togoLFAf_41_clean)
### 

togoLFAf42 <- as.data.frame(cleanlist[[8]]) 

togoLFAf42 %>%
  verify(nrow(togoLFAf42) == 945) %>%
  write_excel_csv(files$togoLFAf_42_clean)

###

togoLFAg41 <- as.data.frame(cleanlist[[9]]) 

togoLFAg41 %>%
  verify(nrow(togoLFAg41) == 1507) %>%
  write_excel_csv(files$togoLFAg_41_clean)

###
togoLFAg42 <- as.data.frame(cleanlist[[10]]) 

togoLFAg42 %>%
  verify(nrow(togoLFAg42) == 1305) %>%
  write_excel_csv(files$togoLFAg_42_clean)

###

togoLFAl41 <- as.data.frame(cleanlist[[11]]) 

togoLFAl41 %>%
  verify(nrow(togoLFAl41) == 1509) %>%
  write_excel_csv(files$togoLFAl_41_clean) 

###
togoLFAl42 <- as.data.frame(cleanlist[[12]]) 

togoLFAl42 %>%
  verify(nrow(togoLFAl42) == 1187) %>%
  write_excel_csv(files$togoLFAl_42_clean) 

###
togoMBAc41 <- as.data.frame(cleanlist[[13]]) 

togoMBAc41 %>%
  verify(nrow(togoMBAc41) == 1513) %>%
  write_excel_csv(files$togoMBAc_41_clean)

###
togoMBAc42 <- as.data.frame(cleanlist[[14]]) 

togoMBAc42 %>%
  verify(nrow(togoMBAc42) == 1397) %>%
  write_excel_csv(files$togoMBAc_42_clean)

###
togoMBAp41 <- as.data.frame(cleanlist[[15]]) 

togoMBAp41 %>%
  verify(nrow(togoMBAp41) == 1513) %>%
  write_excel_csv(files$togoMBAp_41_clean)

###
togoMBAp42 <- as.data.frame(cleanlist[[16]]) 

togoMBAp42 %>%
  verify(nrow(togoMBAp42) == 1397) %>%
  write_excel_csv(files$togoMBAp_42_clean)

# done