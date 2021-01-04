#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JH
# Maintainer(s): JH
# License: 2021, GPL v3 or later
#
# -----------------------------------------------------------
# multicountrytrachseroproject/model/src/model.R

#####################################################################
## Modeling code for                                               ##
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
## modified from                                                   ##
## CODE FOR FITTING MODELS TO CROSS-SECTIONAL ANTIBODY TITRE DATA  ##
##                                                                 ##
## Michael White                                                   ##
## Institut Pasteur                                                ##
## michael.white@pasteur.fr                                        ##
##                                                                 ##
#####################################################################

# load libraries
pacman::p_load(
  "here", "MASS", "compiler",
  "binom", "coda", "readr",
  "janitor", "purrr", "tidyverse", "assertr"
)

# creates a list of all filepaths to data to be read into the script
files <- list(
  drc1_Ct694_cleanmod = here("model/input/drc1_Ct694_cleanmod.csv"),
  drc1_LFA_cleanmod = here("model/input/drc1_LFA_cleanmod.csv"),
  drc1_MBA_cleanmod = here("model/input/drc1_MBA_cleanmod.csv"),
  drc2_Ct694_cleanmod = here("model/input/drc2_Ct694_cleanmod.csv"),
  drc2_LFA_cleanmod = here("model/input/drc2_LFA_cleanmod.csv"),
  drc2_MBA_cleanmod = here("model/input/drc2_MBA_cleanmod.csv"),
  togoLFAf_41_cleanmod = here("model/input/togoLFAf_41_cleanmod.csv"),
  togoLFAf_42_cleanmod = here("model/input/togoLFAf_42_cleanmod.csv"),
  togoLFAg_41_cleanmod = here("model/input/togoLFAg_41_cleanmod.csv"),
  togoLFAg_42_cleanmod = here("model/input/togoLFAg_42_cleanmod.csv"),
  togoLFAl_41_cleanmod = here("model/input/togoLFAl_41_cleanmod.csv"),
  togoLFAl_42_cleanmod = here("model/input/togoLFAl_42_cleanmod.csv"),
  togoMBAc_41_cleanmod = here("model/input/togoMBAc_41_cleanmod.csv"),
  togoMBAc_42_cleanmod = here("model/input/togoMBAc_42_cleanmod.csv"),
  togoMBAp_41_cleanmod = here("model/input/togoMBAp_41_cleanmod.csv"),
  togoMBAp_42_cleanmod = here("model/input/togoMBAp_42_cleanmod.csv")
)

stopifnot(length(files) == 16)

## creates a list of all files as connections
fileslist <- list(
  files$drc1_Ct694_cleanmod, files$drc1_LFA_cleanmod,
  files$drc1_MBA_cleanmod, files$drc2_Ct694_cleanmod,
  files$drc2_LFA_cleanmod, files$drc2_MBA_cleanmod,
  files$togoLFAf_41_cleanmod, files$togoLFAf_42_cleanmod,
  files$togoLFAg_41_cleanmod, files$togoLFAg_42_cleanmod,
  files$togoLFAl_41_cleanmod, files$togoLFAl_42_cleanmod,
  files$togoMBAc_41_cleanmod, files$togoMBAc_42_cleanmod,
  files$togoMBAp_41_cleanmod, files$togoMBAp_42_cleanmod
)

stopifnot(length(fileslist) == 16)

# creates a list called dfs, containing all 16 dataframes created from csvs
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

names(dfs) <- df_names

## using dfs, we extract each df,run it through the model,
## and save and export result to the plot task ##

## 1.1 MODEL
# start i loop
for (i in seq_along(dfs)) {

  # messages for the user to keep them aware of model progress
  start_time <- Sys.time()
  print(paste0("Modelling for dataset ", names(dfs)[i], " has now begun..."))

  # set seed for reproducibility of results
  set.seed(22315)
  seed <- 22315

  df <- as.data.frame(pluck(dfs, i))

  # define the model_M1 function
  model_M1 <- function(a, par_M1) {
    lambda_0 <- par_M1[1]
    rho <- par_M1[2]

    SP_prop <- (lambda_0 / (lambda_0 + rho)) * (1 - exp(-(lambda_0 + rho) * a))

    SP_prop
  } # end model_M1 function

  model_M1 <- cmpfun(model_M1, options = list(optimize = 3))

  ###################################################
  ## 1.2 LIKELIHOOD

  # define the loglike_M1 function
  loglike_M1 <- function(par_M1) {
    SP_model <- sapply(df[, 1], model_M1, par_M1 = par_M1)

    loglike <- df[, 3] * log(SP_model) +
      (1 - df[, 3]) * log(1 - SP_model)

    sum(loglike)
  } # end the loglike_M1 function

  loglike_M1 <- cmpfun(loglike_M1, options = list(optimize = 3))

  ###################################################
  ## 1.3 PRIOR

  LARGE <- 1e10 ## Large value for rejecting parameters with prior

  # define the prior_M1 function
  prior_M1 <- function(par) {
    lambda_0 <- par[1]
    rho <- par[2]

    ###################################################
    ## Uniform prior on lambda_0 ~ U(0,100)

    # start lambda_0 loop
    if (lambda_0 > 0 && lambda_0 < 100) {
      prior_lambda_0 <- log(1 / 100)
    } else {
      prior_lambda_0 <- -LARGE
    } # end lamba_0 loop

    ###################################################
    ## Uniform prior on rho ~ U(0,10)

    # start rho loop
    if (rho > 0 && rho < 10) {
      prior_rho <- log(1 / 10)
    } else {
      prior_rho <- -LARGE
    } # end rho loop

    prior <- prior_lambda_0 + prior_rho

    prior
  } # end prior_M1 loop

  # define the prior_M1 function

  prior_M1 <- cmpfun(prior_M1, options = list(optimize = 3))

  ### MCMC

  N_mcmc <- 10000 ## Number of MCMC iterations

  MCMC_accept <- 0 ## Track the MCMC acceptance rate

  ###################################################
  ## 2.2 Prepare object for MCMC fitting

  N_par <- 2 ## Number of parameters

  MCMC_par <- matrix(NA, nrow = N_mcmc, ncol = N_par + 2)
  colnames(MCMC_par) <- c("lambda_0", "rho", "loglike", "prior")

  ###################################################
  ## 2.3 Implement MCMC iterations

  par_MC <- c(0.1, 0.1) ## (lambda_0, rho)
  par_MCp1 <- rep(NA, N_par) ## Parameter update vector

  Sigma_MC <- c(0.01, 0.01) ## Standard deviation of proposal distribution

  prior_MC <- prior_M1(par_MC)

  loglike_MC <- loglike_M1(par_MC) + prior_MC

  # start mc loop
  for (mc in 1:N_mcmc)
  {
    ###################################################
    ## Propose new parameter

    par_MCp1[1] <- rnorm(n = 1, mean = par_MC[1], sd = Sigma_MC[1])
    par_MCp1[2] <- rnorm(n = 1, mean = par_MC[2], sd = Sigma_MC[2])


    ###################################################
    ## Only proceed if not rejected by the prior

    prior_MCp1 <- prior_M1(par_MCp1)

    # start prior_Mcp1 loop
    if (prior_MCp1 > -0.5 * LARGE) {
      ###################################################
      ## Calculate log-likelihood and use Metropolis-Hastings
      ## algorithm to accept or reject

      loglike_MCp1 <- loglike_M1(par_MCp1) + prior_MCp1

      log_prob <- min(loglike_MCp1 - loglike_MC, 0)

      # start logprob loop
      if (log(runif(1)) < log_prob) {
        par_MC <- par_MCp1

        loglike_MC <- loglike_MCp1
        prior_MC <- prior_MCp1

        MCMC_accept <- MCMC_accept + 1
      } # end logprob loop
    } # end prior_MCp1 loop

    ###################################################

    MCMC_par[mc, 1:N_par] <- par_MC
    MCMC_par[mc, N_par + 1] <- loglike_MC
    MCMC_par[mc, N_par + 2] <- prior_MC
  } # end mc loop

  ###################################################
  ## 3.1 Caluclate the median of the posteriors
  ##     and the model prediction for these
  ##     parameters

  MCMC_burn <- MCMC_par[floor(0.2 * nrow(MCMC_par)):(nrow(MCMC_par) - 1), ]

  age_seq <- seq(from = 0, to = 70, by = 0.2)

  par_median <- apply(X = MCMC_burn[, 1:N_par], MARGIN = 2, FUN = median)

  M1_predict_median <- sapply(age_seq, model_M1, par = par_median)

  ###################################################
  ## 3.2 Take N_sam equally spaced samples from
  ##     the posterior and caluclate the median
  ##     posterior prediction

  N_sam <- 500
  sam_seq <- round(seq(from = 1, to = nrow(MCMC_burn), length = N_sam))


  M1_predict <- matrix(NA, nrow = N_sam, ncol = length(age_seq))
  # begin l loop
  for (l in 1:N_sam)
  {
    M1_predict[l, ] <- sapply(age_seq, model_M1, par = MCMC_burn[sam_seq[l], 1:N_par])
  } # end l loop

  M1_quant <- matrix(NA, nrow = 3, ncol = length(age_seq))

  # begin m loop
  for (m in 1:length(age_seq))
  {
    M1_quant[, m] <- quantile(M1_predict[, m], prob = c(0.025, 0.5, 0.975), seed = seed)
  } # end m loop

  # obtain and export SCR estimates for each site (med, lcl, ucl)

  SCR <- as.data.frame(t(quantile(MCMC_burn[, 1], prob = c(0.5, 0.025, 0.975)))) %>%
    transmute(
      lamba = as.numeric(`50%`),
      ucl_est = as.numeric(`2.5%`),
      lcl_est = as.numeric(`97.5%`)
    )

  write_excel_csv(SCR,
    quote = FALSE,
    path =
      here(paste("model/output/", names(dfs)[i], "_scr_ests.csv",
        sep = ""
      ))
  )

  ###################################################
  # join the model estimates (med, lcl, ucl) and save as new dataframe

  M1_quant_df <- as.data.frame(t(M1_quant)) %>%
    transmute(
      medest = as.numeric(V2),
      high95_est = as.numeric(V3),
      low95_est = as.numeric(V1)
    )

  M1_df <- M1_quant_df %>%
    mutate(rownum = as.numeric(rownames(M1_quant_df)))

  age_seq_df <- as.data.frame(age_seq)
  age_seq_df <- age_seq_df %>%
    mutate(rownum = as.numeric(rownames(M1_quant_df)))

  # export data just for the 1-9 year olds
  model_ests_df <- as.data.frame(left_join(M1_df,
    age_seq_df,
    by = "rownum"
  )) %>%
    filter(age_seq <= 9.0)

  # iterate over each df containing model estimates,
  # name it appropriately, and export it to the plot task } # end i loop

  write_excel_csv(model_ests_df,
    quote = FALSE,
    file =
      here(paste("plot/input/", names(dfs)[i], "_model_ests.csv", sep = ""))
  )

  # message to let the user know that each iteration has completed
  print(paste0("Modelling for dataset ", names(dfs)[i], " has completed successfully."))
} # end i loop

# mark the end of all modeling, let the user know how long it took overall
end_time <- Sys.time()
print(paste0("All modelling complete at ", end_time))
total_time <- end_time - start_time
print(paste0(
  "Running all models took ",
  round(total_time, 2), " minutes to run to completion"
))

# done
