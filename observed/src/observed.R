
#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License: 2020, EICC, GPL v2 or later
#
# -----------------------------------------------------------
# multicountrytrachseroproject/clean/src/clean.R

## Code for obtaining seroprevalence proportions by age group and  ##
## binomial confidence intervals for                               ##
## "Comparison of multiple tests for determination of              ##
## seroconversion rates to the Chlamydia trachomatis antigen       ##
## Pgp3: a multi-country analysis"                                 ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##   
## but please maintain appropriate accreditation)                  ##
##                                                                 ##   
## Jessica Randall                                                 ##
## github.com/jessmrandall                                         ##

## this script reads in cleaned dataframes in each sheet to 
## calculate prevalence proportions for each age group and 95% 
## binomial confidence estimates sheets containing the prevalences 
## and CI estimates are exported to the plot task's input folder 
## for plotting age seroprevalence curves


# load libraries 
pacman::p_load("here","readr", "janitor", "tidyverse", 
               "assertr", "purrr", "binom")

# specify file structure of inputs and outputs

files <- list(
  
  drc1_Ct694_clean = here("/observed/input/DRC1Ct694_clean.csv"),
  drc1_LFA_clean = here("/observed/input/DRC1LFA_clean.csv"),
  drc1_MBA_clean = here("/observed/input/DRC1MBA_clean.csv"),
  drc2_Ct694_clean = here("/observed/input/DRC1Ct694_clean.csv"),
  drc2_LFA_clean = here("/observed/input/DRC1LFA_clean.csv"),
  drc2_MBA_clean = here("/observed/input/DRC1MBA_clean.csv"),
  togoLFAf_41_clean = here("/observed/input/TogoLFAfield_40001_clean.csv"),
  togoLFAf_42_clean = here("/observed/input/TogoLFAfield_40002_clean.csv"),
  togoLFAg_41_clean = here("/observed/input/TogoLFAgold_40001_clean.csv"),
  togoLFAg_42_clean = here("/observed/input/TogoLFAgold_40002_clean.csv"), 
  togoLFAl_41_clean = here("/observed/input/TogoLFAlatex_40001_clean.csv"), 
  togoLFAl_42_clean = here("/observed/input/TogoLFAlatex_40002_clean.csv"),
  togoMBAc_41_clean = here("/observed/input/TogoMBAct694_40001_clean.csv"),
  togoMBAc_42_clean = here("/observed/input/TogoMBAct694_40002_clean.csv"),
  togoMBAp_41_clean = here("/observed/input/TogoMBAct694_40001_clean.csv"),
  togoMBAp_42_clean = here("/observed/input/TogoMBAct694_40002_clean.csv"),
  
  drc1_ct694_obs_plot = here("/plot/input/DRC1Ct694_obs_plot_df.csv"),
  drc1_LFA_obs_plot = here("/plot/input/DRC1LFA_obs_plot_df.csv"),
  drc1_MBA_obs_plot = here("/plot/input/DRC1MBA_obs_plot_df.csv"),
  drc2_ct694_obs_plot = here("/plot/input/DRC2Ct694_obs_plot_df.csv"),
  drc2_LFA_obs_plot = here("/plot/input/DRC2Ct694_obs_plot_df.csv"),
  drc2_MBA_obs_plot = here("/plot/input/DRC2MBA_obs_plot_df.csv"),
  togoLFAf_41_obs_plot = here("/plot/input/togoLFAf_40001_obs_plot_df.csv"),
  togoLFAf_42_obs_plot = here("/plot/input/togoLFAf_40002_obs_plot_df.csv"),
  togoLFAg_41_obs_plot = here("/plot/input/togoLFAgold_40001_obs_plot_df.csv"),
  togoLFAg_42_obs_plot = here("/plot/input/togoLFAgold_40002_obs_plot_df.csv"),
  togoLFAl_41_obs_plot = here("/plot/input/togoLFAlatex_40001_obs_plot_df.csv"),
  togoLFAl_42_obs_plot = here("/plot/input/togoLFAlatex_40002_obs_plot_df.csv"),
  togoMBAc_41_obs_plot = here("/plot/input/togoMBACt694_40001_obs_plot_df.csv"),
  togoMBAc_42_obs_plot = here("/plot/input/togoMBACt694_40004_obs_plot_df.csv"),
  togoMBAp_41_obs_plot = here("/plot/input/togoMBAPgp3_40001_obs_plot_df.csv"),
  togoMBAp_41_obs_plot = here("/plot/input/togoMBAPgp3_40002_obs_plot_df.csv")
)

stopifnot(is_empty(files) != TRUE & length(files) == 32)

#set seed for reproducibility
seed <- 22315
set.seed(22315)

########################################################################
# reads in clean data and checks the number of rows for accuracy
# could also be a loop but i'm not sure how to individually name the 
# resulting dataframes yet

### DRC ###
###)
drc1_ct694 <- as.data.frame(read_csv(files[[1]]))
drc1_ct694 <- drc1_ct694  %>%
  verify(nrow(drc1_ct694) == 1496)
### 
drc1_LFA <- as.data.frame(read_csv(files[[2]]))
drc1_LFA <- drc1_LFA %>%
  verify(nrow(drc1_LFA) == 1494)
###
drc1_MBA <- as.data.frame(read_csv(files[[3]])) 
drc1_MBA <- drc1_MBA %>%
  verify(nrow(drc1_MBA) == 1496)
###
drc2_Ct694 <- as.data.frame(read_csv(files[[4]]))
drc2_Ct694  <- drc2_Ct694 %>%
  verify(nrow(drc2_Ct694) == 1496)
###
drc2_LFA <- as.data.frame(read_csv(files[[5]])) 
drc2_LFA <- drc2_LFA %>%
  verify(nrow(drc2_LFA) == 1494)
###
drc2_MBA <- as.data.frame(read_csv(files[[6]])) 
drc2_MBA <- drc2_MBA %>%
  verify(nrow(drc2_MBA) == 1496)
###

###Togo###
###
togoLFAf41 <- as.data.frame(read_csv(files[[7]])) 
togoLFAf41 <- togoLFAf41 %>%
  verify(nrow(togoLFAf41) == 972)

### 
togoLFAf42 <- as.data.frame(read_csv(files[[8]])) 
togoLFAf42 <- togoLFAf42 %>%
  verify(nrow(togoLFAf42) == 945)
###
togoLFAg41 <- as.data.frame(read_csv(files[[9]])) 
togoLFAg41 <- togoLFAg41 %>%
  verify(nrow(togoLFAg41) == 1507)
###
togoLFAg42 <- as.data.frame(read_csv(files[[10]])) 
togoLFAg42 <- togoLFAg42 %>%
  verify(nrow(togoLFAg42) == 1305)
###
togoLFAl41 <- as.data.frame(read_csv(files[[11]])) 
togoLFAl41 <- togoLFAl41 %>%
  verify(nrow(togoLFAl41) == 1509)
###
togoLFAl42 <- as.data.frame(read_csv(files[[12]])) 
togoLFAl42 <- togoLFAl42 %>%
  verify(nrow(togoLFAl42) == 1187)
###
togoMBAc41 <- as.data.frame(read_csv(files[[13]])) 
togoMBAc41<- togoMBAc41 %>%
  verify(nrow(togoMBAc41) == 1513)
###
togoMBAc42 <- as.data.frame(read_csv(files[[14]])) 
togoMBAc42 <- togoMBAc42 %>%
  verify(nrow(togoMBAc42) == 1397)
###
togoMBAp41 <- as.data.frame(read_csv(files[[15]])) 
togoMBAp41 <- togoMBAp41 %>%
  verify(nrow(togoMBAp41) == 1513)
###
togoMBAp42 <- as.data.frame(read_csv(files[[16]])) 
togoMBAp42 <- togoMBAp42 %>%
  verify(nrow(togoMBAp42) == 1397) 
#############################################

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
#############################################################
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
############################################################

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
  mutate(age = as.factor(row.names(drc1_MBA_sp_bins)))

#min(drc1_MBA_sp_bins$med)
#max(drc1_MBA_sp_bins$med)
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
######################################################

###DRC###
## ct694
drc2_ct694_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(drc2_Ct694$age > age_bins[i] & 
                   drc2_Ct694$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc2_Ct694[index,3]
  
  drc2_ct694_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}

#check that no data are missing 
stopifnot(not_na(drc2_ct694_sp_bins) == TRUE) 

# check the minimum and maximum values in the first column are reproducible
# create a new column called age from the row numbers to join with age bin data
drc2_ct694_sp_bins <-drc2_ct694_sp_bins %>%
  assert(within_bounds(0.08823529, 0.4067797), med) %>%
  mutate(age = as.factor(row.names(drc2_ct694_sp_bins)))

# create country and test specific age bin data frame
drc2_ct694_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(drc2_ct694_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
drc2_ct694_obs<- left_join(drc2_ct694_sp_bins, drc2_ct694_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(drc2_ct694_obs, files$drc2_ct694_obs_plot)
###################################################################
###
# drc2_LFA 
drc2_LFA_sp_bins <- data.frame(med = numeric(0), 
                               low_95 = numeric(0), 
                               high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(drc2_LFA$age > age_bins[i] & 
                   drc2_LFA$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc2_LFA[index,3]
  
  drc2_LFA_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                             length(temp_AB), 
                                                             method="wilson", 
                                                             seed = seed)
                                               [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(drc2_LFA_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
drc2_LFA_sp_bins <-drc2_LFA_sp_bins %>%
  assert(within_bounds(0.0760183,0.3704397),low_95) %>%
  mutate(age = as.factor(row.names(drc2_LFA_sp_bins)))

# create country and test specific age bin data frame
drc2_LFA_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(drc2_LFA_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
drc2_LFA_obs<- left_join(drc2_LFA_sp_bins, drc2_LFA_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(drc2_LFA_obs, files$drc2_LFA_obs_plot)
####################################################
# drc2_MBA 
drc2_MBA_sp_bins <- data.frame(med = numeric(0), 
                               low_95 = numeric(0), 
                               high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(drc2_MBA$age > age_bins[i] & 
                   drc2_MBA$age <= age_bins[i+1] ) 
  
  temp_AB  <- drc2_MBA[index,3]
  
  drc2_MBA_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                             length(temp_AB), 
                                                             method="wilson", 
                                                             seed = seed)
                                               [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(drc2_MBA_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
drc2_MBA_sp_bins <-drc2_MBA_sp_bins %>%
  assert(within_bounds(0.1944542,0.4884959), high_95) %>%
  mutate(age = as.factor(row.names(drc2_MBA_sp_bins)))

# create country and test specific age bin data frame
drc2_MBA_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(drc2_MBA_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
drc2_MBA_obs<- left_join(drc2_MBA_sp_bins, drc2_MBA_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(drc2_MBA_obs, files$drc2_MBA_obs_plot)
######################################################################

### Togo ###
# LFAf41
togoLFAf41_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoLFAf41$age > age_bins[i] & 
                   togoLFAf41$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoLFAf41[index,3]
  
  togoLFAf41_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}

#check that no data are missing 
stopifnot(not_na(togoLFAf41_sp_bins) == TRUE) 

# check the minimum and maximum values in the first column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoLFAf41_sp_bins <-togoLFAf41_sp_bins %>%
  mutate(age = as.factor(row.names(togoLFAf41_sp_bins)))

# create country and test specific age bin data frame
togoLFAf41_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoLFAf41_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoLFAf41_obs<- left_join(togoLFAf41_sp_bins, togoLFAf41_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoLFAf41_obs, files$togoLFAf_41_obs_plot)
#############################################################
###
# Togo LFA field 40002
togoLFAf42_sp_bins <- data.frame(med = numeric(0), 
                               low_95 = numeric(0), 
                               high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoLFAf42$age > age_bins[i] & 
                   togoLFAf42$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoLFAf42[index,3]
  
  togoLFAf42_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                             length(temp_AB), 
                                                             method="wilson", 
                                                             seed = seed)
                                               [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoLFAf42_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoLFAf42_sp_bins <-togoLFAf42_sp_bins %>%
  assert(within_bounds(0,0.01582264),low_95) %>%
  mutate(age = as.factor(row.names(togoLFAf42_sp_bins)))

# create country and test specific age bin data frame
togoLFAf42_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoLFAf42_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoLFAf42_obs<- left_join(togoLFAf42_sp_bins, togoLFAf42_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoLFAf42_obs, files$togoLFAf_42_obs_plot)
############################################################
###
# Togo LFA gold 40001
togoLFAg41_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoLFAg41$age > age_bins[i] & 
                   togoLFAg41$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoLFAg41[index,3]
  
  togoLFAg41_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoLFAg41_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoLFAg41_sp_bins <-togoLFAg41_sp_bins %>%
  mutate(age = as.factor(row.names(togoLFAg41_sp_bins)))

# create country and test specific age bin data frame
togoLFAg41_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoLFAg41_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoLFAg41_obs<- left_join(togoLFAg41_sp_bins, togoLFAg41_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoLFAg41_obs, files$togoLFAg_41_obs_plot)
######################################################
###
# Togo LFA gold 40002
togoLFAg42_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoLFAg42$age > age_bins[i] & 
                   togoLFAg42$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoLFAg42[index,3]
  
  togoLFAg42_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoLFAg42_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoLFAg42_sp_bins <-togoLFAg42_sp_bins %>%
  assert(within_bounds(0.007352941,0.06451613),med) %>%
  mutate(age = as.factor(row.names(togoLFAg42_sp_bins)))

# create country and test specific age bin data frame
togoLFAg42_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoLFAg42_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoLFAg42_obs<- left_join(togoLFAg42_sp_bins, togoLFAg42_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoLFAg42_obs, files$togoLFAg_42_obs_plot)

####################################################################
# LFAl41
###
# Togo LFA latex 40002
togoLFAl41_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoLFAl41$age > age_bins[i] & 
                   togoLFAl41$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoLFAl41[index,3]
  
  togoLFAl41_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoLFAl41_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoLFAl41_sp_bins <-togoLFAl41_sp_bins %>%
  mutate(age = as.factor(row.names(togoLFAl41_sp_bins)))

# create country and test specific age bin data frame
togoLFAl41_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoLFAl41_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoLFAl41_obs<- left_join(togoLFAl41_sp_bins, togoLFAl41_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoLFAl41_obs, files$togoLFAl_41_obs_plot)

######################################################################
# LFAl42
togoLFAl42_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoLFAl42$age > age_bins[i] & 
                   togoLFAl42$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoLFAl42[index,3]
  
  togoLFAl42_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoLFAl42_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoLFAl42_sp_bins <-togoLFAl42_sp_bins %>%
  assert(within_bounds(0,0.07348666),low_95) %>%
  mutate(age = as.factor(row.names(togoLFAl42_sp_bins)))

# create country and test specific age bin data frame
togoLFAl42_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoLFAl42_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoLFAl42_obs<- left_join(togoLFAl42_sp_bins, togoLFAl42_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoLFAl42_obs, files$togoLFAl_42_obs_plot)

###########################################################################
# MBAc41
togoMBAc41_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoMBAc41$age > age_bins[i] & 
                   togoMBAc41$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoMBAc41[index,3]
  
  togoMBAc41_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoMBAc41_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoMBAc41_sp_bins <-togoMBAc41_sp_bins %>%
  assert(within_bounds(0,0.0323647),low_95) %>%
  mutate(age = as.factor(row.names(togoMBAc41_sp_bins)))

# create country and test specific age bin data frame
togoMBAc41_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoMBAc41_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoMBAc41_obs<- left_join(togoMBAc41_sp_bins, togoMBAc41_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoMBAc41_obs, files$togoMBAc_41_obs_plot)

########################################################################
# MBAc42
togoMBAc42_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoMBAc42$age > age_bins[i] & 
                   togoMBAc42$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoMBAc42[index,3]
  
  togoMBAc42_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoMBAc42_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoMBAc42_sp_bins <-togoMBAc42_sp_bins %>%
  mutate(age = as.factor(row.names(togoMBAc42_sp_bins)))

# create country and test specific age bin data frame
togoMBAc42_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoMBAc42_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoMBAc42_obs<- left_join(togoMBAc42_sp_bins, togoMBAc42_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoMBAc42_obs, files$togoMBAc_42_obs_plot)

##############################################################################

# MBAp41
togoMBAp41_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoMBAp41$age > age_bins[i] & 
                   togoMBAp41$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoMBAp41[index,3]
  
  togoMBAp41_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoMBAp41_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoMBAp41_sp_bins <-togoMBAp41_sp_bins %>%
  assert(within_bounds(0,0.05607477),med) %>%
  mutate(age = as.factor(row.names(togoMBAp41_sp_bins)))

# create country and test specific age bin data frame
togoMBAp41_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoMBAp41_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoMBAp41_obs<- left_join(togoMBAp41_sp_bins, togoMBAp41_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoMBAp41_obs, files$togoMBAp_41_obs_plot)

##########################################################################

# MBAp42
togoMBAp42_sp_bins <- data.frame(med = numeric(0), 
                                 low_95 = numeric(0), 
                                 high_95 = numeric(0))

# loop thorugh data to populate the sp_bins into a 9x3 matrix
# containing observed age seroprevalence proportions for each age group and
# 95% confidence intervals
for(i in 1:N_bins)
{
  index <- which(togoMBAp42$age > age_bins[i] & 
                   togoMBAp42$age <= age_bins[i+1] ) 
  
  temp_AB  <- togoMBAp42[index,3]
  
  togoMBAp42_sp_bins[i,] <- as.numeric(as.vector(binom.confint(sum(temp_AB), 
                                                               length(temp_AB), 
                                                               method="wilson", 
                                                               seed = seed)
                                                 [1,4:6]))
}
#check that no data are missing 
stopifnot(not_na(togoMBAp42_sp_bins) == TRUE) 

# check the minimum and maximum values in the second column are reproducible
# create a new column called age from the row numbers to join with age bin data
togoMBAp42_sp_bins <-togoMBAp42_sp_bins %>%
  mutate(age = as.factor(row.names(togoMBAp42_sp_bins)))

# create country and test specific age bin data frame
togoMBAp42_age_bins <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(togoMBAp42_sp_bins$age))

#merge observed prevalence proprtions, confidence intervals, age bins and age 
# for plotting 
togoMBAp42_obs<- left_join(togoMBAp42_sp_bins, togoMBAp42_age_bins, by = "age") %>%
  mutate(age = as.numeric(age))

#export data
write_excel_csv(togoMBAp42_obs, files$togoMBAp_42_obs_plot)

# done