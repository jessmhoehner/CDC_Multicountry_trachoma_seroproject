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
pacman::p_load("MASS", "compiler", "binom", "coda", "readr", "janitor")


# specify inputs and outputs



## 0.1 Read in data

AB_data <- read.csv( file="TogoLFAfield_40002.csv", header=TRUE)


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
 
loglike_M1 <- function( par_M1 )
{
	SP_model <- sapply( AB_data[,1], model_M1, par_M1=par_M1)

	loglike <- AB_data[,3]*log(SP_model) + (1-AB_data[,3])*log(1-SP_model)

      sum( loglike )
}

loglike_M1 <- cmpfun(loglike_M1, options=list(optimize=3))


###################################################
## 1.3 PRIOR

LARGE = 1e10     ## Large value for rejecting parameters with prior
 
prior_M1 <- function( par )
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

N_mcmc       <- 10000      ## Number of MCMC iterations

MCMC_accept <- 0           ## Track the MCMC acceptance rate

#################################################
## 2.2 Prepare object for MCMC fitting

N_par <- 2      ## Number of parameters 
 
MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=N_par+2)
colnames(MCMC_par) <- c("lambda_0", "rho", "loglike", "prior")

 
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
	## export as matrices

	MCMC_par[mc,1:N_par] <- par_MC
	MCMC_par[mc,N_par+1] <- loglike_MC
	MCMC_par[mc,N_par+2] <- prior_MC
}



#########################################################
## 2.4 Examine MCMC chains

# move to plot task and plot in tidyverse
 
par(mfrow=c(3,3))

for(j in 1:N_par)
{
	#####################################
	## PANEL j: MCMC chain

	plot(x=1:N_mcmc, y=MCMC_par[,j], 
	pch=19, col="grey", cex=0.25,
	xlab="MCMC iteration", ylab=colnames(MCMC_par)[j], 
	main=paste("MCMC chain: ", colnames(MCMC_par)[j], " (ESS = ", 
	           round(effectiveSize(MCMC_par[,j]),0), ")",  sep="") )
}

#####################################
## PANEL 3: MCMC log-likelihood

plot(x=1:N_mcmc, y=MCMC_par[,3], 
pch=19, col="grey", cex=0.25,
ylim=quantile(MCMC_par[,3], prob=c(0.01,1)),
xlab="MCMC iteration", ylab="log-likelihood", 
main="log-likelihood" )



#########################################################
## 2.5 Examine posterior distribution
 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]

for(j in 1:N_par)
{
	#####################################
	## PANEL j: MCMC posterior

	DEN <- density( MCMC_burn[,j] )
	
	QUANT <- quantile( MCMC_burn[,j], prob=c(0.025, 0.5, 0.975) )

	plot(x=DEN$x, y=DEN$y, type='l',
	xlim=c(0, max(DEN$x)),
	xlab=colnames(MCMC_par)[j], ylab="", 
	main=paste("Posterior profile: ", colnames(MCMC_par)[j], sep="") )

	
	low_index  = which(DEN$x<QUANT[1])
	mid_index  = intersect( which(DEN$x>=QUANT[1]), which(DEN$x<=QUANT[3]) )
	high_index = which(DEN$x>QUANT[3])

	polygon( x=c( DEN$x[low_index], rev(DEN$x[low_index]) ),
		   y=c( rep(0,length(low_index)), rev(DEN$y[low_index]) ), 
	             col="pink")

	polygon( x=c( DEN$x[mid_index], rev(DEN$x[mid_index]) ),
		   y=c( rep(0,length(mid_index)), rev(DEN$y[mid_index]) ), 
	         col="grey")

	polygon( x=c( DEN$x[high_index], rev(DEN$x[high_index]) ),
		   y=c( rep(0,length(high_index)), rev(DEN$y[high_index]) ), 
		   col="pink")

	points(x=rep(QUANT[2],2), y=c(0,max(DEN$y)), type='l', lty="dashed", lwd=2)
}


plot.new()


#########################################################
## 2.6 Auto-correlation

for(j in 1:N_par)
{
	autocorr.plot( MCMC_par[,j], auto.layout=FALSE, 
                     main=paste("Auto-correlation: ", colnames(MCMC_par)[j]) )
}

plot.new()



paste( "Acceptance rate = ", round(100*MCMC_accept/N_mcmc,2), "%", sep="" )

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
	M1_quant[,j] = quantile( M1_predict[,j], prob=c(0.025, 0.5, 0.975) )
}

## move to plot task ## 
###############################################
## 3.3 Plot data and model prediction
 
par(mfrow=c(1,1))


plot(x=age_bins_mid, y=SP_bins[,1], 
pch=15, cex=2,
xlim=c(0,10), ylim=c(0,1),
xlab="", ylab="", 
main=""  )


for(i in 1:N_bins)
{
	arrows(x0=age_bins_mid[i], y0=SP_bins[i,2], 
             x1=age_bins_mid[i], y1=SP_bins[i,3], 
             length=0.03, angle=90, code=3, col="black", lwd=1)	
}



points(x=age_seq, y=M1_quant[2,], 
type='l', lwd=3, col="blue")

polygon(x=c(age_seq, rev(age_seq)), 
y=c( M1_quant[1,], rev(M1_quant[3,]) ),
col=rgb(0/256,0/256,256/256,0.2), border=NA)



points(x=age_seq, y=M1_predict_median, 
type='l', lwd=3, col="blue", lty="dashed")

quantile( MCMC_burn[,1], prob=c(0.5, 0.025, 0.975) )