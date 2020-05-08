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
  togoMBAp_41_obs_plot = here("/plot/input/togoMBAPgp3_40002_obs_plot_df.csv"),
  
  drc1_ct694_modelests = here("plot/input/drc1_Ct694_model_ests_df.csv"), 
  drc1_LFA_modelests = here("plot/input/drc1_LFA_model_ests_df.csv"),
  drc1_MBA_modelests = here("plot/input/drc1_MBA_model_ests_df.csv"),
  drc2_ct694_modelests = here("plot/input/drc2_Ct694_model_ests_df.csv"), 
  drc2_LFA_modelests = here("plot/input/drc2_LFA_model_ests_df.csv"),
  drc2_MBA_modelests = here("plot/input/drc2_MBA_model_ests_df.csv"),
  togoLFAf_41_modelests = here("/plot/input/togoLFAf_40001_model_ests_df.csv"),
  togoLFAf_42_modelests = here("/plot/input/togoLFAf_40002_model_ests_df.csv"),
  togoLFAg_41_modelests = here("/plot/input/togoLFAgold_40001_model_ests_df.csv"),
  togoLFAg_42_modelests = here("/plot/input/togoLFAgold_40002_model_ests_df.csv"),
  togoLFAl_41_modelests = here("/plot/input/togoLFAlatex_40001_model_ests_df.csv"),
  togoLFAl_42_modelests = here("/plot/input/togoLFAlatex_40002_model_ests_df.csv"),
  togoMBAc_41_modelests = here("/plot/input/togoMBACt694_40001_model_ests_df.csv"),
  togoMBAc_42_modelests = here("/plot/input/togoMBACt694_40004_model_ests_df.csv"),
  togoMBAp_41_modelests = here("/plot/input/togoMBAPgp3_40001_model_ests_df.csv"),
  togoMBAp_41_modelests = here("/plot/input/togoMBAPgp3_40002_model_ests_df.csv"),
  
  drc1_ct694_plot = here("plot/input/drc1_Ct694_plot.png"), 
  drc1_LFA_plot = here("plot/input/drc1_LFA_plot.png"),
  drc1_MBA_plot = here("plot/input/drc1_MBA_plot.png"),
  drc2_ct694_plot = here("plot/input/drc2_Ct694_plot.png"), 
  drc2_LFA_plot = here("plot/input/drc2_LFA_plot.png"),
  drc2_MBA_plot = here("plot/input/drc2_MBA_plot.png"),
  togoLFAf_41_plot = here("/plot/input/togoLFAf_40001_plot.png"),
  togoLFAf_42_plot = here("/plot/input/togoLFAf_40002_plot.png"),
  togoLFAg_41_plot = here("/plot/input/togoLFAgold_40001_plot.png"),
  togoLFAg_42_plot = here("/plot/input/togoLFAgold_40002_plot.png"),
  togoLFAl_41_plot = here("/plot/input/togoLFAlatex_40001_plot.png"),
  togoLFAl_42_plot = here("/plot/input/togoLFAlatex_40002_plot.png"),
  togoMBAc_41_plot = here("/plot/input/togoMBACt694_40001_plot.png"),
  togoMBAc_42_plot = here("/plot/input/togoMBACt694_40004_plot.png"),
  togoMBAp_41_plot = here("/plot/input/togoMBAPgp3_40001_plot.png"),
  togoMBAp_41_plot = here("/plot/input/togoMBAPgp3_40002_plot.png"),
  fullfacetplot= here("/plot/output/fullfacteplot.png")
  )

stopifnot(is_empty(files) != TRUE & length(files) == 49)

pd <- position_dodge(0.1) # move them .05 to the left and right

### plot observed and estimated cases based on SCR model

# read in data as list and save as dfs, combine

# creates a list called dfs, containing all dataframes created from csvs
# to be read in 
inputs <- list(
  drc1_ct694_obs_plot = here("plot/input/DRC1Ct694_obs_plot_df.csv"),
  drc1_LFA_obs_plot = here("plot/input/DRC1LFA_obs_plot_df.csv"),
  drc1_MBA_obs_plot = here("plot/input/DRC1MBA_obs_plot_df.csv"),
  drc2_ct694_obs_plot = here("plot/input/DRC2Ct694_obs_plot_df.csv"),
  drc2_LFA_obs_plot = here("plot/input/DRC2Ct694_obs_plot_df.csv"),
  drc2_MBA_obs_plot = here("plot/input/DRC2MBA_obs_plot_df.csv"),
  togoLFAf_41_obs_plot = here("plot/input/togoLFAf_40001_obs_plot_df.csv"),
  togoLFAf_42_obs_plot = here("plot/input/togoLFAf_40002_obs_plot_df.csv"),
  togoLFAg_41_obs_plot = here("plot/input/togoLFAgold_40001_obs_plot_df.csv"),
  togoLFAg_42_obs_plot = here("plot/input/togoLFAgold_40002_obs_plot_df.csv"),
  togoLFAl_41_obs_plot = here("plot/input/togoLFAlatex_40001_obs_plot_df.csv"),
  togoLFAl_42_obs_plot = here("plot/input/togoLFAlatex_40002_obs_plot_df.csv"),
  togoMBAc_41_obs_plot = here("plot/input/togoMBACt694_40001_obs_plot_df.csv"),
  togoMBAc_42_obs_plot = here("plot/input/togoMBACt694_40002_obs_plot_df.csv"),
  togoMBAp_41_obs_plot = here("plot/input/togoMBAPgp3_40001_obs_plot_df.csv"),
  togoMBAp_42_obs_plot = here("plot/input/togoMBAPgp3_40002_obs_plot_df.csv"),
  
  drc1_ct694_modelests = here("plot/input/drc1_Ct694_model_ests_df.csv"), 
  drc1_LFA_modelests = here("plot/input/drc1_LFA_model_ests_df.csv"),
  drc1_MBA_modelests = here("plot/input/drc1_MBA_model_ests_df.csv"),
  drc2_ct694_modelests = here("plot/input/drc2_Ct694_model_ests_df.csv"), 
  drc2_LFA_modelests = here("plot/input/drc2_LFA_model_ests_df.csv"),
  drc2_MBA_modelests = here("plot/input/drc2_MBA_model_ests_df.csv"),
  
  togoLFAf_41_modelests = here("plot/input/togoLFAf_41_model_ests_df.csv"),
  togoLFAf_42_modelests = here("plot/input/togoLFAf_42_model_ests_df.csv"),
  togoLFAg_41_modelests = here("plot/input/togoLFAg_41_model_ests_df"),
  togoLFAg_42_modelests = here("plot/input/togoLFAg_42_model_ests_df"),
  togoLFAl_41_modelests = here("plot/input/togoLFAl_41_model_ests_df"),
  togoLFAl_42_modelests = here("plot/input/togoLFAl_42_model_ests_df"),
  togoMBAc_41_modelests = here("plot/input/togoMBAc_41_model_ests_df"),
  togoMBAc_42_modelests = here("plot/input/togoMBAc_42_model_ests_df"),
  togoMBAp_41_modelests = here("plot/input/togoMBAp_41_model_ests_df"),
  togoMBAp_42_modelests = here("plot/input/togoMBAp_42_model_ests_df")

)
stopifnot(is_empty(inputs) == FALSE & length(inputs) == 32)

test <-list(inputs$drc1_ct694_obs_plot, inputs$drc1_ct694_modelests)

dfs <- lapply(inputs, function(x) {
  
  df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
    clean_names()
  
}

)

dfs[[1:16]]


# reads in the list of inputs which has two other lists, observedprevs and 
# modelests, checks som expected parameters, and creates a list called dfs

dfs <- lapply(inputs, function(x) {
  
  df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
  clean_names()
  
  obs_df <- df[1:16,]
  
  ests_df <- df[17:32,]
  
  obs_df <- obs_df  %>%
    verify(ncol(obs_df) == 5) %>%
    verify(is.na(obs_df) == FALSE) %>%
    transmute(med = med, 
              high95_obs = high_95,
              low95_obs = low_95,
              age = age, 
              age_bins_mid = age_bins_mid)
  
  ests_df <- ests_df  %>%
    verify(ncol(ests_df) == 5) %>%
    verify(is.na(ests_df) == FALSE) %>%
    transmute(medest = medest, 
              high95_est = high95_est,
              low95_est - low95_est,
              rownum = rownum, 
              age_seq = age_seq)
})

# add names for each df in the list corresponding to appropriate names for each
# spreadheet, in this case country, associated unit and assay information, 
# and whether they describe observed prevalence or model estimate data

df_names <- c("drc1_Ct694_obs" , "drc1_LFA_obs", "drc1_MBA_obs", 
"drc2_Ct694_obs", "drc2_LFA_obs", 
              "drc2_MBA_obs", "togoLFAf_41_obs", "togoLFAf_42_obs", 
			  "togoLFAg_41_obs", 
              "togoLFAg_42_obs", "togoLFAl_41_obs", "togoLFAl_42_obs", 
			  "togoMBAc_41_obs", 
              "togoMBAc_42_obs", "togoMBAp_41_obs", "togoMBAp_42_obs",
			  "drc1_Ct694_obs" , "drc1_LFA_obs", "drc1_MBA_obs", 
"drc2_Ct694_modelests", "drc2_LFA_modelests", 
              "drc2_MBA_modelests", "togoLFAf_41_modelests", 
			  "togoLFAf_42_modelests", 
			  "togoLFAg_41_modelests", 
              "togoLFAg_42_modelests", "togoLFAl_41_modelests", 
			  "togoLFAl_42_modelests", 
			  "togoMBAc_41_modelests", 
              "togoMBAc_42_modelests", "togoMBAp_41_modelests", "togoMBAp_42_modelests",
			  
			  )
			 
names(dfs) <- df_names



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
         subtitle = "Manono, DRC (Ct694) with 95% Confidence Interval Error 
         Bars") +
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
