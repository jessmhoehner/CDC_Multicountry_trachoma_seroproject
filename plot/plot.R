## Graphing code for                                               ##
## "Comparison of multiple tests for determination of seroconversion rates to ##
## the Chlamydia trachomatis antigen Pgp3: a multi-country analysis" ##
##                                                                 ##
## Please feel free to share or modify the code as you see fit     ##   
## but please maintain appropriate accreditation)                  ##
##                                                                 ##   
## Jessica Randall                                                 ##
## github.com/jessmrandall                                         ##

# load libraries 
pacman::p_load("here", "readr", "janitor", 
               "tidyverse", "assertr", "ggplot2", "ggthemes", "gridExtra")

# specify file structure of inputs and outputs
# change to the plot/input/ in final runs

files <- list(
  drc1_ct694agespbins = here("/clean/output/DRC1Ct694_agespbin_df.csv"), 
  drc1_ct694agespbins_plot1 = here("/plot/output/DRC1Ct694_ageseroprev.png"),
  drc1_ct694_modelests = 
    here("/plot/input/DRC1Ct694_modelests.csv"))

stopifnot(is_empty(files) != TRUE & length(files) == 3)

pd <- position_dodge(0.1) # move them .05 to the left and right

# DRC1_Ct694 plot 1 
# read in age_sp_bindf
drc1_ct694_plot1_df <- read_csv(files$drc1_ct694agespbins, 
                                col_names = TRUE, na = "NA") %>%
  clean_names()

drc1_ct694_plot1_df <- drc1_ct694_plot1_df %>%
  mutate(rownum = as.numeric(row.names(drc1_ct694_plot1_df)))

drc1_ct694_ests_df <- read_csv(files$drc1_ct694_modelests, 
                                col_names = TRUE, na = "NA") %>%
  clean_names()

model_plot_df <- full_join(drc1_ct694_plot1_df , drc1_ct694_ests_df, 
                           by = "rownum")

#plot age sero prev curve with error bars
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
  labs(title = "Antibody Positivity Across Age Groups", 
       subtitle = "DRC (Manono) Ct694 with 95% Confidence Interval Error Bars") +
  xlab("Age") +
  ylab("Percent Antibody Positive") + 
  ylim(0, 1.0) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)))

ggsave(files$drc1_ct694agespbins_plot1, 
       plot = last_plot(),dpi = 600,limitsize = TRUE)




# drc1_Ct694 
# drc1_LFA 
# drc1_MBA 
# drc2_Ct694 
# drc2_LFA 
# drc2_MBA 


#########################################################
## 2.4 Examine MCMC chains

# move to plot task and plot in tidyverse
# these files come from the model task

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

# done 
