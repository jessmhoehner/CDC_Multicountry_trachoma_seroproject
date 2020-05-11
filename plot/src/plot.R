#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License: 2020, EICC, GPL v2 or later
#
# -----------------------------------------------------------
# # multicountrytrachseroproject/plot/src/plot.R

####################################################################
## Plotting code for                                               ##
## "Comparison of multiple tests for determination of              ##
## seroconversion rates to the Chlamydia trachomatis antigen       ##
## Pgp3: a multi-country analysis"                                 ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##   
## but please maintain appropriate accreditation)                  ##
##                                                                 ##   
## Jessica Randall                                                 ##
## github.com/jessmrandall                                         ##
#####################################################################

# load libraries 
pacman::p_load("here", "readr", "janitor", 
               "tidyverse", "assertr", "ggplot2")

# specify file structure of inputs and outputs
# change to the plot/input/ in final runs

files <- list(
  drc1_CT694_obs = here("/plot/input/drc1_Ct694_obs_df.csv"),
  drc1_LFA_obs = here("/plot/input/drc1_LFA_obs_df.csv"),
  drc1_MBA_obs = here("/plot/input/drc1_MBA_obs_df.csv"),
  drc2_CT694_obs = here("/plot/input/drc2_Ct694_obs_df.csv"),
  drc2_LFA_obs = here("/plot/input/drc2_LFA_obs_df.csv"),
  drc2_MBA_obs = here("/plot/input/drc2_MBA_obs_df.csv"),
  togoLFAf_41_obs = here("/plot/input/togoLFAf_41_obs_df.csv"),
  togoLFAf_42_obs = here("/plot/input/togoLFAf_42_obs_df.csv"),
  togoLFAg_41_obs = here("/plot/input/togoLFAg_41_obs_df.csv"),
  togoLFAg_42_obs = here("/plot/input/togoLFAg_42_obs_df.csv"),
  togoLFAl_41_obs = here("/plot/input/togoLFAl_41_obs_df.csv"),
  togoLFAl_42_obs = here("/plot/input/togoLFAl_42_obs_df.csv"),
  togoMBAc_41_obs = here("/plot/input/togoMBAc_41_obs_df.csv"),
  togoMBAc_42_obs = here("/plot/input/togoMBAc_42_obs_df.csv"),
  togoMBAp_41_obs = here("/plot/input/togoMBAp_41_obs_df.csv"),
  togoMBAp_42_obs = here("/plot/input/togoMBAp_42_obs_df.csv"),
  
  drc1_CT694_mod = here("plot/input/drc1_Ct694_model_ests_df.csv"),
  drc1_LFA_mod = here("plot/input/drc1_LFA_model_ests_df.csv"),
  drc1_MBA_mod = here("plot/input/drc1_MBA_model_ests_df.csv"),
  drc2_CT694_mod = here("plot/input/drc2_Ct694_model_ests_df.csv"),
  drc2_LFA_mod = here("plot/input/drc2_LFA_model_ests_df.csv"),
  drc2_MBA_mod = here("plot/input/drc2_MBA_model_ests_df.csv"),
  togoLFAf_41_mod = here("plot/input/togoLFAf_41_model_ests_df.csv"),
  togoLFAf_42_mod = here("plot/input/togoLFAf_42_model_ests_df.csv"),
  togoLFAg_41_mod = here("plot/input/togoLFAg_41_model_ests_df.csv"),
  togoLFAg_42_mod = here("plot/input/togoLFAg_42_model_ests_df.csv"), 
  togoLFAl_41_mod = here("plot/input/togoLFAl_41_model_ests_df.csv"), 
  togoLFAl_42_mod = here("plot/input/togoLFAl_42_model_ests_df.csv"),
  togoMBAc_41_mod = here("plot/input/togoMBAc_41_model_ests_df.csv"),
  togoMBAc_42_mod = here("plot/input/togoMBAc_42_model_ests_df.csv"),
  togoMBAp_41_mod = here("plot/input/togoMBAp_41_model_ests_df.csv"),
  togoMBAp_42_mod = here("plot/input/togoMBAp_42_model_ests_df.csv")

)

stopifnot(is_empty(files) != TRUE & length(files) == 32)

pd <- position_dodge(0.1) # move them .05 to the left and right

## creates a list of all files as connections
obslist <- list(files$drc1_CT694_obs,files$drc1_LFA_obs,files$drc1_MBA_obs,
                files$drc2_CT694_obs,files$drc2_LFA_obs,files$drc2_MBA_obs,
                files$togoLFAf_41_obs,files$togoLFAf_42_obs, files$togoLFAg_41_obs,
                files$togoLFAg_42_obs,files$togoLFAl_41_obs,files$togoLFAl_42_obs,
                files$togoMBAc_41_obs,files$togoMBAc_42_obs,files$togoMBAp_41_obs,
                files$togoMBAp_42_obs)

stopifnot(length(obslist) == 16)
stopifnot(obslist %in% files == TRUE)

# creates a list called obsdfs, containing all 16 dataframes created from obs data
obsdfs <- lapply(obslist, function(x) {
  
  x_df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
    clean_names()
  
  x_df  %>%
    verify(ncol(x_df) == 5 & nrow(x_df) == 9)  %>%
    verify(is.na(x_df) == FALSE) %>%
    transmute(med = med, 
              low_95 = low_95, 
              high_95 = high_95, 
              age = age, 
              age_bins_mid = age_bins_mid)
  
  x_df <- x_df %>%
    mutate(rownum = as.numeric(row.names(x_df)))
  
  })

# creates a list called moddfs, containing all 16 dataframes created from modeled data
modlist <- list(files$drc1_CT694_mod,files$drc1_LFA_mod,files$drc1_MBA_mod,
                files$drc2_CT694_mod,files$drc2_LFA_mod,files$drc2_MBA_mod,
                files$togoLFAf_41_mod,files$togoLFAf_42_mod, files$togoLFAg_41_mod,
                files$togoLFAg_42_mod,files$togoLFAl_41_mod,files$togoLFAl_42_mod,
                files$togoMBAc_41_mod,files$togoMBAc_42_mod,files$togoMBAp_41_mod,
                files$togoMBAp_42_mod)

stopifnot(length(modlist) == 16)
stopifnot(modlist %in% files == TRUE)

# creates a list called moddfs, containing all 16 dataframes created from mod data
moddfs <- lapply(modlist, function(x) {
  
  x_df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
    clean_names()
  
  x_df  %>%
    verify(ncol(x_df) == 5 & nrow(x_df) ==46) %>%
    verify(is.na(x_df) == FALSE) %>%
    transmute(medest = medest, 
              low95_est = low95_est, 
              high95_est = high95_est, 
              age_seq = age_seq, 
              rownum = rownum)
  
})

# add names for each df in the list corresponding to appropriate names for each
# spreadheet, in this case country and associated unit and assay information
# this is used for subtitles and unique file names

df_names <- c("Manono CT694 MBA", "Manono LFA Latex", "Manono Pgp3 MBA", 
              "Nyemba CT694 MBA", "Nyemba LFA Latex", "Nyemba Pgp3 MBA", 
              "Keran LFA (Field)", "Anie LFA (Field)", "Keran LFA (Gold)", 
          "Anie LFA (Gold)", "Keran LFA (Latex)", "Anie LFA (Latex)", 
          "Keran CT694 MBA", "Anie CT694 MBA", "Keran Pgp3 MBA", 
          "Anie Pgp3 MBA")

names(obsdfs) <- df_names
names(moddfs) <- df_names

## using dfs, we extract each df and combine them as needed for each graph
# initialize empty dfs
df_obs <- data.frame(matrix(0, 9, 6))
df_mod <- data.frame(matrix(0, 46, 5))
plot_df <- data.frame(matrix(0, 46, 10))

# make age seroprevalence graphs using the observed data and modeled fit estimates

#start j loop
suppressWarnings(
  for (j in seq_along(obslist)){
  
  #messages for the user to keep them aware of model progress
  print(paste0("Creating plotting dataset for ",names(obsdfs)[j]))
    
  df_obs <- as.data.frame(pluck(obsdfs, j))
  df_mod <- as.data.frame(pluck(moddfs, j))
  plot_df <- full_join(df_obs , df_mod, by = "rownum")
  
  #plot age sero prev curves, observed and estimated with error bars and 95%CI
  age_seroprev_plot<-ggplot(plot_df, aes(age, med, age_bins_mid)) +
      geom_pointrange(aes(ymin=low_95, 
                          ymax=high_95)) +
      geom_errorbar(aes(ymin=low_95,
                        ymax=high_95), 
                    colour="black", 
                    width=.1, 
                    position=pd) +
      geom_line(aes(age_seq, medest), 
                color = "blue", 
                position=pd) +
      geom_ribbon(aes(age_seq,
                      ymin = low95_est, 
                      ymax = high95_est), 
                  fill = "blue", 
                  position=pd, 
                  alpha = 0.2) + 
      theme_classic() +
      labs(title = "Proportion of Antibody Positivity Across Age Groups", 
           subtitle = 
             paste("plot of ",names(obsdfs)[j],
                   "with 95% Confidence Interval Error Bars", sep = "")) +
      xlab("Age") +
      ylab("Proportion seropositive") + 
      ylim(0, 1.0) +
      scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)) +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x  = element_text(size=14),
          axis.text.y  = element_text(size=14))
  
  # save each graph individually
  ggsave(filename = here(paste("plot/output/",names(obsdfs)[j],
                               "_ageseroprev.png", sep = "")), 
         plot = last_plot(),
         device = "png",
         dpi = 600)
  
  #message to let the user know that each iteration has completed
  print(paste0("Plot created for ",names(obsdfs)[j]))
  
}

)


# Suppressed warnings because 
# "Removed 37 rows containing missing values (geom_pointrange)" is because
# there are empty rows in the observed data where the modeled data uses 
# the age_seq variable, I'm hoping to figure out a more elegant solution but this
# plots the correct data for now

# done 
