# TREX1_Reaction_Simulator.R

# A script to simulate TREX1-like reactions over time

##########################################################################

## Environment

rm(list=ls())
setwd('~/path/to/file/')

## Define Parameters

# General
savefile='save.file.name'

# Simulation
time=60*60 # s, total simulation time
save=1 # s, save file time-step
dt=1e-3 # s, integration time-step

# Reaction Components
beta=100 # [nt]/[S], excisable nucleotides per polynucleotide
E.0=7.5e-9 # M, initial enzyme concentration
M.0=0e-9 # M, initial mutant concentration
S.0=0.83e-9 # M, initial substrate concentration
P.0=0 # M, initial product concentration
ES.0=0 # M, initial enzyme-substrate concentration
EP.0=0 # M, initial enzyme-product concentration
MS.0=0 # M, initial mutant-substrate concentration
MP.0=0 # M, intial mutant-product concentration

# Wild-type Enzyme Constants
psi=0.75 # processivity factor
k.n1=0.050 # s^-1 # [ES]->[E][S] rate constant
k.1=1e9 # s^-1 # [E][S]->[ES] rate constant
k.2=10 # s^-1 # [ES]->[E][P]([S]) rate constant + [ES]->[ES][P] rate constant
k.n3=0 # s^-1 # [EP]->[E][P] rate constant
k.3=0 # s^-1 # [E][P]->[EP] rate constant

# Mutant Enzyme Constants
phi=psi # processivity factor 
k.na=k.n1 # s^-1 # [MS]->[M][S] rate constant
k.a=k.1 # s^-1 # [M][S]->[MS] rate constant
k.b=0 # s^-1 # [MS]->[M][P]([S]) rate constant + [MS]->[MS][P] rate constant
k.nc=k.n3 # s^-1 # [MP]->[M][P] rate constant
k.c=k.3 # s^-1 # [M][P]->[MP] rate constant

#### BEGIN SCRIPT

## Generate Initial Conditions

# Create Vectors
t=0
E=E.0
M=M.0
S=matrix(0,nrow=beta,ncol=(time/save+1))
S[beta,1]=S.0
P=P.0
ES=matrix(0,nrow=beta,ncol=(time/save+1))
ES[beta,1]=ES.0
MS=matrix(0,nrow=beta,ncol=(time/save+1))
MS[beta,1]=MS.0
EP=EP.0
MP=MP.0

# Initial Concentrations
t.now=0
E.now=E.0
M.now=M.0
S.now=rep(0,times=beta)
S.now[beta]=S.0
P.now=P.0
ES.now=rep(0,times=beta)
ES.now[beta]=ES.0
MS.now=rep(0,times=beta)
MS.now[beta]=MS.0
EP.now=EP.0
MP.now=MP.0

# Initial Rates
d.E=k.n3*EP.now+k.n1*sum(ES.now)+k.2*(1-psi)*sum(ES.now[2:beta])+k.2*ES.now[1]-k.3*E.now*P.now-k.1*E.now*sum(S.now)
d.M=k.nc*MP.now+k.na*sum(MS.now)+k.b*(1-phi)*sum(MS.now[2:beta])+k.b*MS.now[1]-k.c*M.now*P.now-k.a*M.now*sum(S.now)
d.P=k.n3*EP.now+k.nc*MP.now+k.2*sum(ES.now)+k.b*sum(MS.now)-k.3*E.now*P.now-k.c*M.now*P.now
d.S=k.n1*ES.now+k.na*MS.now+k.2*(1-psi)*c(ES.now[2:beta],0)+k.b*(1-phi)*c(MS.now[2:beta],0)-k.1*E.now*S.now-k.a*M.now*S.now
d.ES=k.1*E.now*S.now+k.2*psi*c(ES.now[2:beta],0)-(k.n1+k.2)*ES.now
d.EP=k.3*E.now*P.now-k.n3*EP.now
d.MS=k.a*M.now*S.now+k.b*phi*c(MS.now[2:beta],0)-(k.na+k.b)*MS.now
d.MP=k.c*M.now*P.now-k.nc*MP.now

## Start Simulation

for (i in 1:(time/dt)){
  t.now=t.now+dt
  E.now=E.now+d.E*dt
  M.now=M.now+d.M*dt
  S.now=S.now+d.S*dt
  P.now=P.now+d.P*dt
  ES.now=ES.now+d.ES*dt
  EP.now=EP.now+d.EP*dt
  MS.now=MS.now+d.MS*dt
  MP.now=MP.now+d.MP*dt
  d.E=k.n3*EP.now+k.n1*sum(ES.now)+k.2*(1-psi)*sum(ES.now[2:beta])+k.2*ES.now[1]-k.3*E.now*P.now-k.1*E.now*sum(S.now)
  d.M=k.nc*MP.now+k.na*sum(MS.now)+k.b*(1-phi)*sum(MS.now[2:beta])+k.b*MS.now[1]-k.c*M.now*P.now-k.a*M.now*sum(S.now)
  d.P=k.n3*EP.now+k.nc*MP.now+k.2*sum(ES.now)+k.b*sum(MS.now)-k.3*E.now*P.now-k.c*M.now*P.now
  d.S=k.n1*ES.now+k.na*MS.now+k.2*(1-psi)*c(ES.now[2:beta],0)+k.b*(1-phi)*c(MS.now[2:beta],0)-k.1*E.now*S.now-k.a*M.now*S.now
  d.ES=k.1*E.now*S.now+k.2*psi*c(ES.now[2:beta],0)-(k.n1+k.2)*ES.now
  d.EP=k.3*E.now*P.now-k.n3*EP.now
  d.MS=k.a*M.now*S.now+k.b*phi*c(MS.now[2:beta],0)-(k.na+k.b)*MS.now
  d.MP=k.c*M.now*P.now-k.nc*MP.now
  
  if (sum(seq(0,time/dt,save/dt)==i)>=1){
    t=append(t,t.now,after = length(t))
    E=append(E,E.now,after = length(E))
    M=append(M,M.now,after = length(M))
    S[,(t.now/save+1)]=S.now
    P=append(P,P.now,after = length(P))
    ES[,(t.now/save+1)]=ES.now
    EP=append(EP,EP.now,after = length(EP))
    MS[,(t.now/save+1)]=MS.now
    MP=append(MP,MP.now,after = length(MP))
    show(paste("Progress = ",round(i*dt/time*100,1),"%",sep=""))
  }
}

#show((P+EP+MP)[round(t)==20*60])
#show(sum(rev((S+ES+MS)[,1200])))
#show(rev((S+ES+MS)[,1200]))

## Save Simulation Data

save(list=ls(),file=paste(savefile,"_Results.RData",sep = ""))

#### END SCRIPT

## Graph Results

# Relative Signal
plot(t/60,(1-((P+EP+MP)/beta/S.0)),type='l',col='black',ylab = 'Relative Signal',xlab = 't (min)',main='Signal vs Time',ylim=c(0,1))

# Product Levels
plot(t/60,(P+EP+MP),type='l',col='black',ylab = 'Excised dNMPs (M)',xlab = 't (min)',main='Product vs Time',ylim = c(0,6e-8))



