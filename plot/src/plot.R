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
  drc1_ct694_observed = here("/plot/input/DRC1Ct694_observed_df.csv"), 
  drc1_ct694_modelests = 
    here("/plot/input/DRC1Ct694_modelests.csv"),
  
  drc1_ct694agespbins_plot1 = here("/plot/output/DRC1Ct694_ageseroprev.png"))

stopifnot(is_empty(files) != TRUE & length(files) == 3)

pd <- position_dodge(0.1) # move them .05 to the left and right

### plot observed and estimated cases based on SCR model

# DRC1_Ct694 plot 1 
# read in observed data
drc1_ct694_plot1_df <- read_csv(files$drc1_ct694_observed, 
                                col_names = TRUE, na = "NA") %>%
  clean_names()

#prepare for joining with ests
drc1_ct694_plot1_df <- drc1_ct694_plot1_df %>%
  mutate(rownum = as.numeric(row.names(drc1_ct694_plot1_df)))

# read in model estimates
drc1_ct694_ests_df <- read_csv(files$drc1_ct694_modelests, 
                               col_names = TRUE, na = "NA") %>%
  clean_names()

# join observed and estimated cases
model_plot_df <- full_join(drc1_ct694_plot1_df , drc1_ct694_ests_df, 
                           by = "rownum")

#plot age sero prev curves, observed and estimated with error bars and 95%CI
(drc_1_ct694_1 <- ggplot(model_plot_df, aes(age, med, age_bins_mid)) +
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
         subtitle = "Manono, DRC (Ct694) with 95% Confidence Interval Error Bars") +
    xlab("Age") +
    ylab("Proportion seropositive") + 
    ylim(0, 1.0) +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size=14),
        axis.text.y  = element_text(size=14))

ggsave(files$drc1_ct694agespbins_plot1, 
       plot = last_plot(),dpi = 600,limitsize = TRUE)
