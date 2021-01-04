#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JH
# Maintainer(s): JH
# License: 2021, GPL v3 or later
#
# -----------------------------------------------------------
# multicountrytrachseroproject/clean/src/clean.R

#####################################################################
## Code for obtaining seroprevalence proportions by age group and  ##
## binomial confidence intervals for                               ##
## "Comparison of platforms for testing antibodies to Chlamydia    ##
## trachomatis antigens: data from the Democratic Republic of the  ##
## Congo and Togo"                                                 ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##
## but please maintain appropriate accreditation)                  ##
##                                                                 ##
## Jessica Hoehner                                                 ##
## github.com/jessmhoehner                                         ##
#####################################################################

## this script reads in cleaned dataframes in each sheet to
## calculate prevalence proportions for each age group and 95%
## binomial confidence estimates sheets containing the prevalences
## and CI estimates are exported to the plot task's input folder
## for plotting age seroprevalence curves


# load libraries
pacman::p_load("here", "readr", "janitor", "tidyverse", "assertr", "binom")

# specify file structure of inputs and outputs

files <- list(
  drc1_Ct694_clean = here("observed/input/drc1_Ct694_cleanobs.csv"),
  drc1_LFA_clean = here("observed/input/drc1_LFA_cleanobs.csv"),
  drc1_MBA_clean = here("observed/input/drc1_MBA_cleanobs.csv"),
  drc2_Ct694_clean = here("observed/input/drc2_Ct694_cleanobs.csv"),
  drc2_LFA_clean = here("observed/input/drc2_LFA_cleanobs.csv"),
  drc2_MBA_clean = here("observed/input/drc2_MBA_cleanobs.csv"),
  togoLFAf_41_clean = here("observed/input/togoLFAf_41_cleanobs.csv"),
  togoLFAf_42_clean = here("observed/input/togoLFAf_42_cleanobs.csv"),
  togoLFAg_41_clean = here("observed/input/togoLFAg_41_cleanobs.csv"),
  togoLFAg_42_clean = here("observed/input/togoLFAg_42_cleanobs.csv"),
  togoLFAl_41_clean = here("observed/input/togoLFAl_41_cleanobs.csv"),
  togoLFAl_42_clean = here("observed/input/togoLFAl_42_cleanobs.csv"),
  togoMBAc_41_clean = here("observed/input/togoMBAc_41_cleanobs.csv"),
  togoMBAc_42_clean = here("observed/input/togoMBAc_42_cleanobs.csv"),
  togoMBAp_41_clean = here("observed/input/togoMBAp_41_cleanobs.csv"),
  togoMBAp_42_clean = here("observed/input/togoMBAp_42_cleanobs.csv")
)

stopifnot(is_empty(files) != TRUE & length(files) == 16)

########################################################################
# reads in clean data and checks the number of rows for accuracy
# could also be a loop but i'm not sure how to individually name the
# resulting dataframes yet
# creates a list called dfs, containing all 16 dataframes created from csvs

fileslist <- list(
  files$drc1_Ct694_clean, files$drc1_LFA_clean,
  files$drc1_MBA_clean, files$drc2_Ct694_clean,
  files$drc2_LFA_clean, files$drc2_MBA_clean,
  files$togoLFAf_41_clean, files$togoLFAf_42_clean,
  files$togoLFAg_41_clean, files$togoLFAg_42_clean,
  files$togoLFAl_41_clean, files$togoLFAl_42_clean,
  files$togoMBAc_41_clean, files$togoMBAc_42_clean,
  files$togoMBAp_41_clean, files$togoMBAp_42_clean
)

stopifnot(length(fileslist) == 16)

# testfiles <- fileslist[1:2]

# must be a list of connections
dfs <- lapply(fileslist, function(x) {
  x_df <- as.data.frame(read_csv(x, col_names = TRUE, na = "NA")) %>%
    clean_names()

  x_df %>%
    verify(ncol(x_df) == 3) %>%
    verify(is.na(x_df) == FALSE) %>%
    transmute(
      age = age,
      titre = titre,
      sero_pos = sero_pos
    )
})

# add names for each df in the list corresponding to appropriate names for each
# spreadheet, in this case country and associated unit and assay information

df_names <- c(
  "drc1_Ct694", "drc1_LFA", "drc1_MBA", "drc2_Ct694", "drc2_LFA",
  "drc2_MBA", "togoLFAf_41", "togoLFAf_42", "togoLFAg_41",
  "togoLFAg_42", "togoLFAl_41", "togoLFAl_42", "togoMBAc_41",
  "togoMBAc_42", "togoMBAp_41", "togoMBAp_42"
)

# df_names_test <- df_names[1:2]

names(dfs) <- df_names

#############################################

# loop though 16 datasets to calculate observed age seroprevalence
# start i loop
for (i in seq_along(dfs)) {
  # set seed for reproducibility of results
  set.seed(22315)
  seed <- 22315

  # messages for the user to keep them aware of model progress
  start_time <- Sys.time()
  print(paste0("Age seroprevalence calculation for dataset ", names(dfs)[i], " has now begun..."))

  df <- as.data.frame(pluck(dfs, i))

  # age bin the data for plotting
  age_bins <- seq(from = 0, to = 9, by = 1)
  age_bins_mid <- seq(from = 0.5, to = 8.5, by = 1)

  N_bins <- length(age_bins) - 1

  # initialize empty dataframe to fill with binomial confidence intervals
  # from observed data
  sp_bins <- data.frame(
    med = numeric(0),
    low_95 = numeric(0),
    high_95 = numeric(0)
  )

  # loop thorugh data to populate the sp_bins into a 9x3 matrix
  # containing observed age seroprevalence proportions for each age group and
  # 95% confidence intervals
  # start k loop
  for (k in 1:N_bins) {
    index <- which(df[, 1] > age_bins[k] & df[, 1] <= age_bins[k + 1])

    temp_AB <- df[index, 3]

    sp_bins[k, ] <- as.numeric(as.vector(binom.confint(sum(temp_AB),
      length(temp_AB),
      method = "wilson",
      seed = seed
    )[1, 4:6]))
  } # close k loop

  # create a new column called age from the row numbers to join with
  # age bin data

  sp_bins <- as.data.frame(sp_bins) %>%
    mutate(age = as.numeric(row.names(sp_bins)))

  # check that no data are missing
  stopifnot(not_na(sp_bins) == TRUE)

  # create country and test specific age bin data frame
  age_bins <- as.data.frame(age_bins_mid) %>%
    filter(age_bins_mid != 9.5) %>%
    mutate(age = as.numeric(sp_bins$age))

  # merge observed prevalence proprtions, confidence intervals, age bins and age
  # for plotting
  obs <- left_join(sp_bins, age_bins, by = "age") %>%
    mutate(age = as.numeric(sp_bins$age))

  # export each df to plot task
  write_excel_csv(obs,
    quote = FALSE, path =
      here(paste("plot/input/", names(dfs)[i], "_obs.csv", sep = ""))
  )

  write_excel_csv(obs,
    quote = FALSE, path =
      here(paste("adjust/input/", names(dfs)[i], "_obs.csv", sep = ""))
  )

  # message to let the user know that each iteration has completed
  print(paste0("Age seroprevalence for dataset ", names(dfs)[i], " has completed successfully."))
} # close i loop

# done