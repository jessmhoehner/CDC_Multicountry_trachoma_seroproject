#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING MODELS TO CROSS-SECTIONAL ANTIBODY TITRE DATA  ##
##                                                                 ##
## Please feel free to share modify the code as you see fit        ##   
## (but please maintain appropriate accreditation)                 ##
##                                                                 ##   
## Michael White                                                   ##
## Institut Pasteur                                                ##
## michael.white@pasteur.fr                                        ##
##                                                                 ##
#####################################################################
#####################################################################

# load libraries 
pacman::p_load("here", "MASS", "compiler", 
               "binom", "coda", "readr", 
               "janitor", "purrr", "tidyverse", "assertr")


# specify inputs and outputs
files <- list(drc1_Ct694_clean = here("/model/input/DRC1Ct694_clean.csv"),
              drc1_ct694agespbins = 
                here("/model/input/DRC1Ct694_agespbin_df.csv"), 
              
              drc1_ct694_modelests = 
              here("/plot/input/DRC1Ct694_modelests.csv"))

stopifnot(is_empty(files) != TRUE & length(files) == 5)

# set random seed
set.seed(22315)            

## Read in data

# DRC

# DRC1_CT694

drc1_ct694_df <- read_csv(files$drc1_Ct694_clean, col_names = TRUE, 
                          na = "NA") %>%
  clean_names()

################################################### 
## 1.1 MODEL   

model_M1 = function(a, par_M1)
{
	lambda_0 <- par_M1[1]
	rho      <- par_M1[2]

	SP_prop <- ( lambda_0/(lambda_0+rho) )*( 1 - exp(-(lambda_0+rho)*a) ) 

	SP_prop
}
 
model_M1 <- cmpfun(model_M1, options=list(optimize=3)) 


###################################################
## 1.2 LIKELIHOOD
## changed where the model refers to column number to specific column name
 
loglike_M1 <- function(par_M1)
{
	SP_model <- sapply( drc1_ct694_df$age, model_M1, par_M1=par_M1)

	loglike <- drc1_ct694_df$sero_pos*log(SP_model) + 
	  (1-drc1_ct694_df$sero_pos)*log(1-SP_model)

      sum(loglike)
}

loglike_M1 <- cmpfun(loglike_M1, options=list(optimize=3))

###################################################
## 1.3 PRIOR

LARGE = 1e10     ## Large value for rejecting parameters with prior
 
prior_M1 <- function(par)
{
	lambda_0 <- par[1]
      rho      <- par[2]

	######################################
	## Uniform prior on lambda_0 ~ U(0,100)

	if( lambda_0>0 && lambda_0<100 )
	{
		prior_lambda_0 <- log(1/100)
	}else{
		prior_lambda_0 <- -LARGE
	}

	######################################
	## Uniform prior on rho ~ U(0,10)

	if( rho>0 && rho<10 )
	{
		prior_rho <- log(1/10)
	}else{
		prior_rho <- -LARGE
	}

	prior <- prior_lambda_0 + prior_rho

	prior
}

prior_M1 <- cmpfun(prior_M1, options=list(optimize=3))

### MCMC

N_mcmc <- 10000            ## Number of MCMC iterations

MCMC_accept <- 0           ## Track the MCMC acceptance rate

#################################################
## 2.2 Prepare object for MCMC fitting

N_par <- 2      ## Number of parameters 
 
MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=N_par+2)
colnames(MCMC_par) <- c("lambda_0", "rho", "loglike", "prior")
# update this later

#########################################################
## 2.3 Implement MCMC iterations

par_MC   <- c(0.1, 0.1)       ## (lambda_0, rho)
par_MCp1 <- rep(NA, N_par)    ## Parameter update vector

Sigma_MC <- c(0.01, 0.01)     ## Standard deviation of proposal distribution

prior_MC   <- prior_M1( par_MC )

loglike_MC <- loglike_M1( par_MC ) + prior_MC

for(mc in 1:N_mcmc)
{
	#######################################
	## Propose new parameter

	par_MCp1[1] <- rnorm(n=1, mean=par_MC[1], sd=Sigma_MC[1])
	par_MCp1[2] <- rnorm(n=1, mean=par_MC[2], sd=Sigma_MC[2])


 	##############################################
	## Only proceed if not rejected by the prior
 
	prior_MCp1 <- prior_M1(par_MCp1)

	if( prior_MCp1 > -0.5*LARGE  )
	{
	 	##############################################
		## Calculate log-likelihood and use Metropolis-Hastings
		## algorithm to accept or reject
 
		loglike_MCp1 <- loglike_M1( par_MCp1 ) + prior_MCp1

		log_prob <- min( loglike_MCp1-loglike_MC, 0 )           
                   
		if( log(runif(1)) < log_prob ) 
		{
			par_MC <- par_MCp1

			loglike_MC  <- loglike_MCp1
			prior_MC    <- prior_MCp1

			MCMC_accept <- MCMC_accept + 1  
		}
	}

	#######################################

	MCMC_par[mc,1:N_par] <- par_MC
	MCMC_par[mc,N_par+1] <- loglike_MC
	MCMC_par[mc,N_par+2] <- prior_MC
}

## keep in model task ##

#############################################
#############################################
##          ##                             ##
##   ####   ##  ###### #####  ###  ######  ##
##  ##  ##  ##    ##   ##    ##      ##    ##
##     ##   ##    ##   ####   ###    ##    ##
##  ##  ##  ##    ##   ##       ##   ##    ##
##   ####   ##    ##   #####  ###    ##    ##
##          ##                             ##
#############################################
#############################################

#############################################
## 3.1 Caluclate the median of the posteriors 
##     and the model prediction for these
##     parameters

MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]

age_seq <- seq(from=0, to=70, by=0.2)

par_median <- apply(X=MCMC_burn[,1:N_par], MARGIN=2, FUN=median)

M1_predict_median <- sapply(age_seq, model_M1, par=par_median)

#############################################
## 3.2 Take N_sam equally spaced samples from 
##     the posterior and caluclate the median
##     posterior prediction

N_sam = 500
sam_seq = round(seq(from=1, to=nrow(MCMC_burn), length=N_sam))


M1_predict = matrix(NA, nrow=N_sam, ncol=length(age_seq))
for(k in 1:N_sam)
{
	M1_predict[k,] = sapply(age_seq, model_M1, par=MCMC_burn[sam_seq[k],1:N_par])
}


M1_quant = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
	M1_quant[,j] = quantile(M1_predict[,j], prob=c(0.025, 0.5, 0.975),seed = seed)
}

# export M1_quant to plot task 
# join all of the model estimated params
# by rownumber and export to plot task

M1_quant_df <- as.data.frame(t(M1_quant)) %>%
  transmute(medest = as.numeric(V2), 
            high95_est = as.numeric(V3), 
            low95_est = as.numeric(V1))
M1_quant_df <- M1_quant_df %>%
  mutate(rownum = as.numeric(rownames(M1_quant_df)))
  
M1_predmed_df <- as.data.frame(M1_predict_median)
M1_predmed_df <- M1_predmed_df %>%
  mutate(rownum = as.numeric(rownames(M1_quant_df)))

M1_df <- left_join(M1_quant_df,M1_predmed_df, by = "rownum")

age_seq_df <- as.data.frame(age_seq)
age_seq_df  <- age_seq_df %>%
  mutate(rownum = as.numeric(rownames(M1_quant_df)))

model_ests_df <- left_join(M1_df,
                           age_seq_df, 
                           by = "rownum") %>%
  filter(age_seq <= 1.6) %>%
  mutate(age = rownum)

stopifnot(nrow(model_ests_df) == 9 & ncol(model_ests_df) == 7)
stopifnot(is_empty(model_ests_df) == FALSE)

write_excel_csv(model_ests_df, files$drc1_ct694_modelests)


## move to plot task ## 
###############################################
## 3.3 Plot data and model prediction
 
par(mfrow=c(1,1))


drc1_ct694_age_spins<- read_csv(files$drc1_ct694agespbins, col_names = TRUE, 
                          na = "NA") %>%
  clean_names()

plot(x=drc1_ct694_age_spins$age_bins_mid, y=drc1_ct694_age_spins$med, 
pch=15, cex=2,
xlim=c(0,10), ylim=c(0,1),
xlab="", ylab="", 
main=""  )


for(i in 1:10)
{
	arrows(x0=drc1_ct694_age_spins$age_bins_mid[i], 
	       y0=drc1_ct694_age_spins$low_95[i], 
	       x1=drc1_ct694_age_spins$age_bins_mid[i], 
	       y1=drc1_ct694_age_spins$high_95[i], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}


# here's where the model fit is added
points(x=model_ests_df$age_seq, y=model_ests_df$medest, 
type='l', lwd=3, col="blue")

# fit line and blue shading
polygon(x=c(model_ests_df$age_seq, rev(model_ests_df$age_seq)), 
y=c(model_ests_df$low95_est, rev(model_ests_df$high95_est) ),
col=rgb(0/256,0/256,256/256,0.2), border=NA)

quantile(MCMC_burn[,1], prob=c(0.5, 0.025, 0.975) )
