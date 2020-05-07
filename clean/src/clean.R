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
  
  drc1_Ct694_clean = here("/model/input/DRC1Ct694_clean.csv"),
  drc1_ct694_obs_plot = here("/plot/input/DRC1Ct694_obs_plot_df.csv"),
  
  drc1_LFA_clean = here("/model/input/DRC1LFA_clean.csv"),
  drc1_LFA_obs_plot = here("/plot/input/DRC1LFA_obs_plot_df.csv"),
  
  drc1_MBA_clean = here("/model/input/DRC1MBA_clean.csv"),
  drc1_MBA_obs_plot = here("/plot/input/DRC1MBA_obs_plot_df.csv"),
  
  drc2_Ct694_clean = here("/model/input/DRC1Ct694_clean.csv"),
  drc2_ct694_obs_plot = here("/plot/input/DRC2Ct694_obs_plot_df.csv"),
  
  drc2_LFA_clean = here("/model/input/DRC1LFA_clean.csv"),
  drc2_LFA_obs_plot = here("/plot/input/DRC2Ct694_obs_plot_df.csv"),
  
  drc2_MBA_clean = here("/model/input/DRC1MBA_clean.csv"),
  drc2_MBA_obs_plot = here("/plot/input/DRC2MBA_obs_plot_df.csv"),
  
  togoLFAf_41_clean = here("/model/input/TogoLFAfield_40001_clean.csv"),
  togoLFAf_41_obs_plot = here("/plot/input/togoLFAf_40001_obs_plot_df.csv"),
  
  togoLFAf_42_clean = here("/model/input/TogoLFAfield_40002_clean.csv"),
  togoLFAf_42_obs_plot = here("/plot/input/togoLFAf_40002_obs_plot_df.csv"),
  
  togoLFAg_41_clean = here("/model/input/TogoLFAgold_40001_clean.csv"),
  togoLFAg_41_obs_plot = here("/plot/input/togoLFAgold_40001_obs_plot_df.csv"),
  
  togoLFAg_42_clean = here("/model/input/TogoLFAgold_40002_clean.csv"), 
  togoLFAg_42_obs_plot = here("/plot/input/togoLFAgold_40002_obs_plot_df.csv"),
  
  togoLFAl_41_clean = here("/model/input/TogoLFAlatex_40001_clean.csv"), 
  togoLFAl_41_obs_plot = here("/plot/input/togoLFAlatex_40001_obs_plot_df.csv"),
  
  togoLFAl_42_clean = here("/model/input/TogoLFAlatex_40002_clean.csv"),
  togoLFAl_42_obs_plot = here("/plot/input/togoLFAlatex_40002_obs_plot_df.csv"),
  
  togoMBAc_41_clean = here("/model/input/TogoMBAct694_40001_clean.csv"),
  togoMBAc_41_obs_plot = here("/plot/input/togoMBACt694_40001_obs_plot_df.csv"),
  
  togoMBAc_42_clean = here("/model/input/TogoMBAct694_40002_clean.csv"),
  togoMBAc_42_obs_plot = here("/plot/input/togoMBACt694_40004_obs_plot_df.csv"),
  
  togoMBAp_41_clean = here("/model/input/TogoMBAct694_40001_clean.csv"),
  togoMBAp_41_obs_plot = here("/plot/input/togoMBAPgp3_40001_obs_plot_df.csv"),
  
  togoMBAp_42_clean = here("/model/input/TogoMBAct694_40002_clean.csv"),
  togoMBAp_41_obs_plot = here("/plot/input/togoMBAPgp3_40002_obs_plot_df.csv")
)

stopifnot(is_empty(files) != TRUE & length(files) == 48)
## Read in data

#set seed for reproducibility
seed = set.seed(22315)

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

#check number of rows (vary by data frame) and export to next task as list#

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

drc2_Ct694_df %>%
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
  write_excel_csv(files$togoLFAf41_clean)
### 

togoLFAf42 <- as.data.frame(cleanlist[[8]]) 

togoLFAf42 %>%
  verify(nrow(togoLFAf42) == 945) %>%
  write_excel_csv(files$togoLFAf42_clean)

###

togoLFAg41 <- as.data.frame(cleanlist[[9]]) 

togoLFAg41 %>%
  verify(nrow(togoLFAg41) == 1507) %>%
  write_excel_csv(files$togoLFAg41_clean)

###
togoLFAg42 <- as.data.frame(cleanlist[[10]]) 

togoLFAg42 %>%
  verify(nrow(togoLFAg42) == 1305) %>%
  write_excel_csv(files$togoLFAg42_clean)

###

togoLFAl41 <- as.data.frame(cleanlist[[11]]) 

togoLFAl41 %>%
  verify(nrow(togoLFAl41) == 1509) %>%
  write_excel_csv(files$togoLFAl41_clean) 

###
togoLFAl42 <- as.data.frame(cleanlist[[12]]) 

togoLFAl42 %>%
  verify(nrow(togoLFAl42) == 1187) %>%
  write_excel_csv(files$togoLFAl42_clean) 

###
togoMBAc41 <- as.data.frame(cleanlist[[13]]) 

togoMBAc41 %>%
  verify(nrow(togoMBAc41) == 1513) %>%
  write_excel_csv(files$togoMBAc41_clean)

###
togoMBAc42 <- as.data.frame(cleanlist[[14]]) 

togoMBAc42 %>%
  verify(nrow(togoMBAc42) == 1397) %>%
  write_excel_csv(files$togoMBAc42_clean)

###
togoMBAp41 <- as.data.frame(cleanlist[[15]]) 

togoMBAp41 %>%
  verify(nrow(togoMBAp41) == 1513) %>%
  write_excel_csv(files$togoMBAp41_clean)

###
togoMBAp42 <- as.data.frame(cleanlist[[16]]) 

togoMBAp42 %>%
  verify(nrow(togoMBAp42) == 1397) %>%
  write_excel_csv(files$togoMBAp42_clean)


###########################################################################

# age bin the data for plotting
age_bins <- seq(from=0, to=9, by=1)
age_bins_mid <- seq(from=0.5, to=8.5, by=1) 

N_bins <- length(age_bins) - 1

#############################################

# initialize empty dataframe to fill with binomial confidence intervals
# from observed data 

###DRC###

## ct694
drc1_ct694_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(drc1_ct694$age > age_bins[i] & 
                   drc1_ct694$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc1_ct694[index,3]
  
  drc1_ct694_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(drc1_ct694_sp_bins) == TRUE) 

# check the minimum and maximum values in the first column are reproducible
# create a new column called age from the row numbers to join with age bin data
drc1_ct694_sp_bins <-drc1_ct694_sp_bins %>%
  assert(within_bounds(0.08823529, 0.4067797), med) %>%
  mutate(age = as.factor(row.names(drc1_ct694_sp_bins)))

# create country and test specific age bin data frame
drc1_ct694_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(drc1_ct694_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
drc1_ct694_obs<- left_join(drc1_ct694_sp_bins, drc1_ct694_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(drc1_ct694_obs, files$drc1_ct694_obs_plot)
###

# drc1_LFA 
drc1_LFA_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(drc1_LFA$age > age_bins[i] & 
                   drc1_LFA$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc1_LFA[index,3]
  
  drc1_LFA_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(drc1_LFA_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
drc1_LFA_sp_bins <-drc1_LFA_sp_bins %>%
  assert(within_bounds(0.0760183,0.3704397),low_95) %>%
  mutate(age = as.factor(row.names(drc1_LFA_sp_bins)))

# create country and test specific age bin data frame
drc1_LFA_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(drc1_LFA_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
drc1_LFA_obs<- left_join(drc1_LFA_sp_bins, drc1_LFA_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(drc1_LFA_obs, files$drc1_LFA_obs_plot)

###

# drc1_MBA 

drc1_MBA_sp_bins <- data.frame(med = numeric(0), 
                               low_95 = numeric(0), 
                               high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(drc1_MBA$age > age_bins[i] & 
                   drc1_MBA$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc1_MBA[index,3]
  
  drc1_MBA_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                             length(temp_AB), 
                                                             method="wilson", 
                                                             seed = seed)
                                               [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(drc1_MBA_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
drc1_MBA_sp_bins <-drc1_MBA_sp_bins %>%
  assert(within_bounds(0.1176471,0.4108108), med) %>%
  mutate(age = as.factor(row.names(drc1_MBA_sp_bins)))

min(drc1_MBA_sp_bins$med)
max(drc1_MBA_sp_bins$med)

# stopped here 5/6/ assert throws an error even though these are the min and max

# create country and test specific age bin data frame
drc1_MBA_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(drc1_MBA_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
drc1_MBA_obs<- left_join(drc1_MBA_sp_bins, drc1_MBA_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(drc1_MBA_obs, files$drc1_MBA_obs_plot)

# drc2_Ct694 
# drc2_LFA 
# drc2_MBA 

# Togo
# LFAf41
# LFAf42
# LFAg41
# LFAg42
# LFAl41
# LFAl42
# MBAc41
# MBAc42
# MBAp41
# MBAp42


# done