# Kcat-Km_Approximation.R

# Monte-Carlo Simulation of hT1 ssDNA Kinetic Variances

##########################################################################################
##########################################################################################

## Environment

rm(list=ls())
#setwd('/path/to/file/')
setwd('~/Desktop/')

## Master Parameters

precision=1e6 # N/A, sample size per simulation
loops=1e0 # N/A, number of simulations
correlation='no' # is sampling data correlated
S.1=15 # nM, low concentration of substrate
S.2=515 # nM, high concentration of substrate
#
mean.1=1.9
sd.1=0.45
mean.2=13
sd.2=1.8

#### BEGIN SCRIPT

### hT1_+/+

Kcat.mean=rep(0,times=loops)
Kcat.sd=rep(0,times=loops)
Km.mean=rep(0,times=loops)
Km.sd=rep(0,times=loops)
CE.mean=rep(0,times=loops)
CE.sd=rep(0,times=loops)

for (i in 1:loops){
  
  ## Define Parameters
  
  CR.1=rgamma(precision,shape=(mean.1/sd.1)^2,rate = (mean.1/sd.1^2))
  CR.2=rgamma(precision,shape=(mean.2/sd.2)^2,rate = (mean.2/sd.2^2))
  if (correlation=='yes'){
    CR.1=CR.1[order(CR.1)]
    CR.2=CR.2[order(CR.2)]
  }
  
  ## Perform Calculations
  
  Kcat = (-S.2 + S.1)*CR.2*CR.1/(S.1*CR.2 - S.2*CR.1)
  Kcat.mean[i] = mean(Kcat)
  Kcat.sd[i]=sd(Kcat)
  Km = (-CR.2 + CR.1)*S.2*S.1/(S.1*CR.2 - S.2*CR.1)
  Km.mean[i] = mean(Km)
  Km.sd[i]=sd(Km)
  CE=Kcat/Km
  CE.mean[i]=mean(CE)
  CE.sd[i]=sd(CE)
}

Kcat.mean=mean(Kcat.mean)
Kcat.sd=median(Kcat.sd)
Km.mean=mean(Km.mean)
Km.sd=median(Km.sd)
CE.mean=mean(CE.mean)
CE.sd=median(CE.sd)

## Save Calculations

Results.WT_WT=list('Kcat.mean'=Kcat.mean,'Kcat.sd'=Kcat.sd,'Km.mean'=Km.mean,'Km.sd'=Km.sd,'CE.mean'=CE.mean,'CE.sd'=CE.sd)
View(Results.WT_WT)
rm(Kcat,Kcat.mean,Kcat.sd,Km,Km.mean,Km.sd,CR.1,CR.2,i,CE,CE.mean,CE.sd)

### Export Data

save(list=ls(),file='Kcat-Km_Results.RData')

#### END SCRIPT

