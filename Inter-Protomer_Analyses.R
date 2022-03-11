# A script to analyze simulations associated with TREX1 Inter-Protomer Communication

#############################################################################
#############################################################################
#############################################################################

##### Environment

rm(list=ls())
setwd('/path/to/files/')

#############################################################################

##### Global RMSD

####################################

### mT1_Apo

# load data

rm(list=ls())
load('mT1-Apo_Input.RData')

# extract relevant 

mT1.Apo.labels=mT1.Apo[['labels']]
mT1.Apo.xyz=mT1.Apo[['xyz']]

rm(mT1.Apo)

# reformat for analysis

mT1.Apo.xyz.trim=mT1.Apo.xyz[,(mT1.Apo.labels[,5]=='protA' | mT1.Apo.labels[,5]=='protB'),,]
rm(mT1.Apo.xyz)

mT1.Apo.xyz.trim.1=matrix(mT1.Apo.xyz.trim[,,,1],ncol = prod(dim(mT1.Apo.xyz.trim)[1:2]),byrow = TRUE)
mT1.Apo.xyz.trim.2=matrix(mT1.Apo.xyz.trim[,,,2],ncol = prod(dim(mT1.Apo.xyz.trim)[1:2]),byrow = TRUE)
mT1.Apo.xyz.trim.3=matrix(mT1.Apo.xyz.trim[,,,3],ncol = prod(dim(mT1.Apo.xyz.trim)[1:2]),byrow = TRUE)
mT1.Apo.xyz.trim.4=matrix(mT1.Apo.xyz.trim[,,,4],ncol = prod(dim(mT1.Apo.xyz.trim)[1:2]),byrow = TRUE)

rm(mT1.Apo.xyz.trim)

# perform analysis

mT1.Apo.rmsd.1=sqrt(rowMeans((mT1.Apo.xyz.trim.1-matrix(mT1.Apo.xyz.trim.1[1,],nrow = nrow(mT1.Apo.xyz.trim.1),ncol = ncol(mT1.Apo.xyz.trim.1),byrow = TRUE))^2))
mT1.Apo.rmsd.2=sqrt(rowMeans((mT1.Apo.xyz.trim.2-matrix(mT1.Apo.xyz.trim.2[1,],nrow = nrow(mT1.Apo.xyz.trim.2),ncol = ncol(mT1.Apo.xyz.trim.2),byrow = TRUE))^2))
mT1.Apo.rmsd.3=sqrt(rowMeans((mT1.Apo.xyz.trim.3-matrix(mT1.Apo.xyz.trim.3[1,],nrow = nrow(mT1.Apo.xyz.trim.3),ncol = ncol(mT1.Apo.xyz.trim.3),byrow = TRUE))^2))
mT1.Apo.rmsd.4=sqrt(rowMeans((mT1.Apo.xyz.trim.4-matrix(mT1.Apo.xyz.trim.4[1,],nrow = nrow(mT1.Apo.xyz.trim.4),ncol = ncol(mT1.Apo.xyz.trim.4),byrow = TRUE))^2))

mT1.Apo.rmsd.all=cbind(mT1.Apo.rmsd.1,mT1.Apo.rmsd.2,mT1.Apo.rmsd.3,mT1.Apo.rmsd.4)*sqrt(3)
mT1.Apo.rmsd.mean=rowMeans(mT1.Apo.rmsd.all)
mT1.Apo.rmsd.sd=apply(mT1.Apo.rmsd.all,MARGIN = 1,sd)

mT1.Apo.rmsd.global=list('all'=mT1.Apo.rmsd.all,'mean'=mT1.Apo.rmsd.mean,'sd'=mT1.Apo.rmsd.sd)
RMSD.Global=list('mT1.Apo'=mT1.Apo.rmsd.global)

# plot check

plot(seq(0,1000,0.2),mT1.Apo.rmsd.mean+mT1.Apo.rmsd.sd,type='l',lty='dotted',xlab = 'Time (ns)',ylab = 'RMSD (angstroms)',main='mT1.Apo',col='grey')
lines(seq(0,1000,0.2),mT1.Apo.rmsd.mean-mT1.Apo.rmsd.sd,lty='dotted',col='grey')
lines(seq(0,1000,0.2),mT1.Apo.rmsd.mean)

# save and cleanup

rm(list=setdiff(ls(),'RMSD.Global'))

####################################

### mT1_4-mer

# load data

load('mT1-4mer_Input.RData')

# extract relevant 

mT1.4mer.labels=mT1.4mer[['labels']]
mT1.4mer.xyz=mT1.4mer[['xyz']]

rm(mT1.4mer)

# reformat for analysis

mT1.4mer.xyz.trim=mT1.4mer.xyz[,(mT1.4mer.labels[,5]=='protA' | mT1.4mer.labels[,5]=='protB'),,]
rm(mT1.4mer.xyz)

mT1.4mer.xyz.trim.1=matrix(mT1.4mer.xyz.trim[,,,1],ncol = prod(dim(mT1.4mer.xyz.trim)[1:2]),byrow = TRUE)
mT1.4mer.xyz.trim.2=matrix(mT1.4mer.xyz.trim[,,,2],ncol = prod(dim(mT1.4mer.xyz.trim)[1:2]),byrow = TRUE)
mT1.4mer.xyz.trim.3=matrix(mT1.4mer.xyz.trim[,,,3],ncol = prod(dim(mT1.4mer.xyz.trim)[1:2]),byrow = TRUE)
mT1.4mer.xyz.trim.4=matrix(mT1.4mer.xyz.trim[,,,4],ncol = prod(dim(mT1.4mer.xyz.trim)[1:2]),byrow = TRUE)

rm(mT1.4mer.xyz.trim)

# perform analysis

mT1.4mer.rmsd.1=sqrt(rowMeans((mT1.4mer.xyz.trim.1-matrix(mT1.4mer.xyz.trim.1[1,],nrow = nrow(mT1.4mer.xyz.trim.1),ncol = ncol(mT1.4mer.xyz.trim.1),byrow = TRUE))^2))
mT1.4mer.rmsd.2=sqrt(rowMeans((mT1.4mer.xyz.trim.2-matrix(mT1.4mer.xyz.trim.2[1,],nrow = nrow(mT1.4mer.xyz.trim.2),ncol = ncol(mT1.4mer.xyz.trim.2),byrow = TRUE))^2))
mT1.4mer.rmsd.3=sqrt(rowMeans((mT1.4mer.xyz.trim.3-matrix(mT1.4mer.xyz.trim.3[1,],nrow = nrow(mT1.4mer.xyz.trim.3),ncol = ncol(mT1.4mer.xyz.trim.3),byrow = TRUE))^2))
mT1.4mer.rmsd.4=sqrt(rowMeans((mT1.4mer.xyz.trim.4-matrix(mT1.4mer.xyz.trim.4[1,],nrow = nrow(mT1.4mer.xyz.trim.4),ncol = ncol(mT1.4mer.xyz.trim.4),byrow = TRUE))^2))

mT1.4mer.rmsd.all=cbind(mT1.4mer.rmsd.1,mT1.4mer.rmsd.2,mT1.4mer.rmsd.3,mT1.4mer.rmsd.4)*sqrt(3)
mT1.4mer.rmsd.mean=rowMeans(mT1.4mer.rmsd.all)
mT1.4mer.rmsd.sd=apply(mT1.4mer.rmsd.all,MARGIN = 1,sd)

mT1.4mer.rmsd.global=list('all'=mT1.4mer.rmsd.all,'mean'=mT1.4mer.rmsd.mean,'sd'=mT1.4mer.rmsd.sd)
RMSD.Global[['mT1.4mer']]=mT1.4mer.rmsd.global

# plot check

plot(seq(0,1000,0.2),mT1.4mer.rmsd.mean+mT1.4mer.rmsd.sd,type='l',lty='dotted',xlab = 'Time (ns)',ylab = 'RMSD (angstroms)',main='mT1.4mer',col='grey')
lines(seq(0,1000,0.2),mT1.4mer.rmsd.mean-mT1.4mer.rmsd.sd,lty='dotted',col='grey')
lines(seq(0,1000,0.2),mT1.4mer.rmsd.mean)

# save and cleanup

rm(list=setdiff(ls(),'RMSD.Global'))

####################################

### mT1_dsDNA

# load data

load('mT1-dsDNA_Input.RData')

# extract relevant 

mT1.dsDNA.labels=mT1.dsDNA[['labels']]
mT1.dsDNA.xyz=mT1.dsDNA[['xyz']]

rm(mT1.dsDNA)

# reformat for analysis

mT1.dsDNA.xyz.trim=mT1.dsDNA.xyz[,(mT1.dsDNA.labels[,5]=='protA' | mT1.dsDNA.labels[,5]=='protB'),,]
rm(mT1.dsDNA.xyz)

mT1.dsDNA.xyz.trim.1=matrix(mT1.dsDNA.xyz.trim[,,,1],ncol = prod(dim(mT1.dsDNA.xyz.trim)[1:2]),byrow = TRUE)
mT1.dsDNA.xyz.trim.2=matrix(mT1.dsDNA.xyz.trim[,,,2],ncol = prod(dim(mT1.dsDNA.xyz.trim)[1:2]),byrow = TRUE)
mT1.dsDNA.xyz.trim.3=matrix(mT1.dsDNA.xyz.trim[,,,3],ncol = prod(dim(mT1.dsDNA.xyz.trim)[1:2]),byrow = TRUE)
mT1.dsDNA.xyz.trim.4=matrix(mT1.dsDNA.xyz.trim[,,,4],ncol = prod(dim(mT1.dsDNA.xyz.trim)[1:2]),byrow = TRUE)

rm(mT1.dsDNA.xyz.trim)

# perform analysis

mT1.dsDNA.rmsd.1=sqrt(rowMeans((mT1.dsDNA.xyz.trim.1-matrix(mT1.dsDNA.xyz.trim.1[1,],nrow = nrow(mT1.dsDNA.xyz.trim.1),ncol = ncol(mT1.dsDNA.xyz.trim.1),byrow = TRUE))^2))
mT1.dsDNA.rmsd.2=sqrt(rowMeans((mT1.dsDNA.xyz.trim.2-matrix(mT1.dsDNA.xyz.trim.2[1,],nrow = nrow(mT1.dsDNA.xyz.trim.2),ncol = ncol(mT1.dsDNA.xyz.trim.2),byrow = TRUE))^2))
mT1.dsDNA.rmsd.3=sqrt(rowMeans((mT1.dsDNA.xyz.trim.3-matrix(mT1.dsDNA.xyz.trim.3[1,],nrow = nrow(mT1.dsDNA.xyz.trim.3),ncol = ncol(mT1.dsDNA.xyz.trim.3),byrow = TRUE))^2))
mT1.dsDNA.rmsd.4=sqrt(rowMeans((mT1.dsDNA.xyz.trim.4-matrix(mT1.dsDNA.xyz.trim.4[1,],nrow = nrow(mT1.dsDNA.xyz.trim.4),ncol = ncol(mT1.dsDNA.xyz.trim.4),byrow = TRUE))^2))

mT1.dsDNA.rmsd.all=cbind(mT1.dsDNA.rmsd.1,mT1.dsDNA.rmsd.2,mT1.dsDNA.rmsd.3,mT1.dsDNA.rmsd.4)*sqrt(3)
mT1.dsDNA.rmsd.mean=rowMeans(mT1.dsDNA.rmsd.all)
mT1.dsDNA.rmsd.sd=apply(mT1.dsDNA.rmsd.all,MARGIN = 1,sd)

mT1.dsDNA.rmsd.global=list('all'=mT1.dsDNA.rmsd.all,'mean'=mT1.dsDNA.rmsd.mean,'sd'=mT1.dsDNA.rmsd.sd)
RMSD.Global[['mT1.dsDNA']]=mT1.dsDNA.rmsd.global

# plot check

plot(seq(0,1000,0.2),mT1.dsDNA.rmsd.mean+mT1.dsDNA.rmsd.sd,type='l',lty='dotted',xlab = 'Time (ns)',ylab = 'RMSD (angstroms)',main='mT1.dsDNA',col='grey')
lines(seq(0,1000,0.2),mT1.dsDNA.rmsd.mean-mT1.dsDNA.rmsd.sd,lty='dotted',col='grey')
lines(seq(0,1000,0.2),mT1.dsDNA.rmsd.mean)

# save and cleanup

rm(list=setdiff(ls(),'RMSD.Global'))

save(RMSD.Global,file = 'RMSD-Global.RData')

rm(list=ls())

#############################################################################

##### AlphaCarbon RMSF

####################################

### mT1_Apo

# load data

rm(list=ls())
load('mT1-Apo_Input.RData')

# extract relevant 

mT1.Apo.labels=mT1.Apo[['labels']]
mT1.Apo.xyz=mT1.Apo[['xyz']]

rm(mT1.Apo)

# reformat for analysis

mT1.Apo.xyz.trim=mT1.Apo.xyz[,((mT1.Apo.labels[,5]=='protB') & mT1.Apo.labels[,2]=='CA'),,]
rm(mT1.Apo.xyz)

# analyze

temp.1=apply(mT1.Apo.xyz.trim,MARGIN = c(1,2,4),FUN = mean)
temp.2=mT1.Apo.xyz.trim
for(i in 1:dim(mT1.Apo.xyz.trim)[3]){
  temp.2[,,i,]=temp.1
}

mT1.Apo.rmsf.AlphaCarbon.all=apply(apply((mT1.Apo.xyz.trim-temp.2)^2,MARGIN = c(2,3,4),FUN = sum),MARGIN = c(1,3),FUN = mean)
row.names(mT1.Apo.rmsf.AlphaCarbon.all)=paste(rep(c('B'),each=231),'_',mT1.Apo.labels[((mT1.Apo.labels[,5]=='protB') & mT1.Apo.labels[,2]=='CA'),4],'-',rep(4:234,times=1),sep='')
mT1.Apo.residues=row.names(mT1.Apo.rmsf.AlphaCarbon.all)
rm(temp.1,temp.2,i)

mT1.Apo.rmsf.AlphaCarbon.mean=rowMeans(mT1.Apo.rmsf.AlphaCarbon.all)
mT1.Apo.rmsf.AlphaCarbon.sd=apply(mT1.Apo.rmsf.AlphaCarbon.all,MARGIN = 1,FUN = sd)

mT1.Apo.rmsf.AlphaCarbon=list('residues'=mT1.Apo.residues,'all'=mT1.Apo.rmsf.AlphaCarbon.all,'mean'=mT1.Apo.rmsf.AlphaCarbon.mean,'sd'=mT1.Apo.rmsf.AlphaCarbon.sd)

RMSF.AlphaCarbon=list('mT1.Apo'=mT1.Apo.rmsf.AlphaCarbon)

# save and cleanup

rm(list=setdiff(ls(),'RMSF.AlphaCarbon'))

####################################

### mT1_dsDNA

# load data

load('mT1-dsDNA_Input.RData')

# extract relevant 

mT1.dsDNA.labels=mT1.dsDNA[['labels']]
mT1.dsDNA.xyz=mT1.dsDNA[['xyz']]

rm(mT1.dsDNA)

# reformat for analysis

mT1.dsDNA.xyz.trim=mT1.dsDNA.xyz[,((mT1.dsDNA.labels[,5]=='protB') & mT1.dsDNA.labels[,2]=='CA'),,]
rm(mT1.dsDNA.xyz)

# analyze

temp.1=apply(mT1.dsDNA.xyz.trim,MARGIN = c(1,2,4),FUN = mean)
temp.2=mT1.dsDNA.xyz.trim
for(i in 1:dim(mT1.dsDNA.xyz.trim)[3]){
  temp.2[,,i,]=temp.1
}

mT1.dsDNA.rmsf.AlphaCarbon.all=apply(apply((mT1.dsDNA.xyz.trim-temp.2)^2,MARGIN = c(2,3,4),FUN = sum),MARGIN = c(1,3),FUN = mean)
row.names(mT1.dsDNA.rmsf.AlphaCarbon.all)=paste(rep(c('B'),each=231),'_',mT1.dsDNA.labels[((mT1.dsDNA.labels[,5]=='protB') & mT1.dsDNA.labels[,2]=='CA'),4],'-',rep(4:234,times=1),sep='')
mT1.dsDNA.residues=row.names(mT1.dsDNA.rmsf.AlphaCarbon.all)
rm(temp.1,temp.2,i)

mT1.dsDNA.rmsf.AlphaCarbon.mean=rowMeans(mT1.dsDNA.rmsf.AlphaCarbon.all)
mT1.dsDNA.rmsf.AlphaCarbon.sd=apply(mT1.dsDNA.rmsf.AlphaCarbon.all,MARGIN = 1,FUN = sd)

mT1.dsDNA.rmsf.AlphaCarbon=list('residues'=mT1.dsDNA.residues,'all'=mT1.dsDNA.rmsf.AlphaCarbon.all,'mean'=mT1.dsDNA.rmsf.AlphaCarbon.mean,'sd'=mT1.dsDNA.rmsf.AlphaCarbon.sd)

RMSF.AlphaCarbon[['mT1.dsDNA']]=mT1.dsDNA.rmsf.AlphaCarbon

# result check

plot(1:231,RMSF.AlphaCarbon[["mT1.Apo"]][["mean"]],type='l',col='black',main='mTREX1 Backbone Fluctuations',ylab = 'RMSF (Å)',xlab='Residues (N → C)',ylim=c(0,50),xaxt='n')
lines(1:231,RMSF.AlphaCarbon[["mT1.dsDNA"]][["mean"]],col='red')
abline(v=c(160.5,172.5),col='green')
abline(v=c(19.5,24.5),col='purple')
abline(v=c(70.5,95.5),col='orange')
legend('topright',legend = c('mT1','mT1 + dsDNA'),fill = c('black','red'),col = c('black','red'))
axis(1,at = 115.5,labels = 'Protomer-B')

# save and cleanup

rm(list=setdiff(ls(),'RMSF.AlphaCarbon'))

save(RMSF.AlphaCarbon,file = 'RMSF-AlphaCarbon.RData')

rm(list = ls())

#############################################################################

##### AlphaCarbon Energy Map

####################################

# load and trim data

rm(list=ls())

load('mT1-Apo_Input.RData')
mT1.Apo.labels=mT1.Apo[['labels']]
mT1.Apo.xyz=mT1.Apo[['xyz']]
mT1.Apo.xyz.alphas=mT1.Apo.xyz[,((mT1.Apo.labels[,5]=='protB') & mT1.Apo.labels[,2]=='CA'),,]
rm(mT1.Apo,mT1.Apo.xyz)

load('mT1-dsDNA_Input.RData')
mT1.dsDNA.labels=mT1.dsDNA[['labels']]
mT1.dsDNA.xyz=mT1.dsDNA[['xyz']]
mT1.dsDNA.xyz.alphas=mT1.dsDNA.xyz[,((mT1.dsDNA.labels[,5]=='protB') & mT1.dsDNA.labels[,2]=='CA'),,]
rm(mT1.dsDNA,mT1.dsDNA.xyz)

# reformat

Alphas.data=matrix(c(mT1.Apo.xyz.alphas,mT1.dsDNA.xyz.alphas),ncol=231*3,byrow = TRUE)
Alphas.key=cbind(rep(c('mT1.Apo','mT1.dsDNA'),each=5001*4),rep(rep(1:4,each=5001),times=2))
colnames(Alphas.key)=c('SIM','RUN')
row.names(Alphas.data)=paste(Alphas.key[,1],'_',Alphas.key[,2],sep='')
EnergyMap.AlphaCarbon.input=list('data'=Alphas.data,'key'=Alphas.key)
rm(list=setdiff(ls(),'EnergyMap.AlphaCarbon.input'))

# analyze

EnergyMap.AlphaCarbon.princomp.raw=princomp(EnergyMap.AlphaCarbon.input[['data']])

EnergyMap.AlphaCarbon.princomp.VarProp=(EnergyMap.AlphaCarbon.princomp.raw[['sdev']]^2)/sum(EnergyMap.AlphaCarbon.princomp.raw[['sdev']]^2)

temp.1=EnergyMap.AlphaCarbon.princomp.raw[['scores']][,1:2]
colnames(temp.1)=c('PC.1','PC.2')
EnergyMap.AlphaCarbon.Maps.mT1.Apo=temp.1[EnergyMap.AlphaCarbon.input[['key']][,1]=='mT1.Apo',]
EnergyMap.AlphaCarbon.Maps.mT1.dsDNA=temp.1[EnergyMap.AlphaCarbon.input[['key']][,1]=='mT1.dsDNA',]
rm(temp.1)

# reconcile

EnergyMap.AlphaCarbon.Maps=list('mT1.Apo'=EnergyMap.AlphaCarbon.Maps.mT1.Apo,'mT1.dsDNA'=EnergyMap.AlphaCarbon.Maps.mT1.dsDNA)
EnergyMap.AlphaCarbon.princomp=list('raw'=EnergyMap.AlphaCarbon.princomp.raw,'VarProp'=EnergyMap.AlphaCarbon.princomp.VarProp)

EnergyMap.AlphaCarbon=list('input'=EnergyMap.AlphaCarbon.input,'princomp'=EnergyMap.AlphaCarbon.princomp,'Maps'=EnergyMap.AlphaCarbon.Maps)
rm(list=setdiff(ls(),'EnergyMap.AlphaCarbon'))

# plot check

library(MASS)

contour(kde2d(EnergyMap.AlphaCarbon[["Maps"]][["mT1.Apo"]][,1],EnergyMap.AlphaCarbon[["Maps"]][["mT1.Apo"]][,2]),nlevels = 10,col = 'black',drawlabels = F,xlim = c(-50,50),ylim = c(-50,50),main='Cα Dynamic Range (Protomer-B)',xlab='PC.1 (25.2%)',ylab='PC.2 (12.4%)')
contour(kde2d(EnergyMap.AlphaCarbon[["Maps"]][["mT1.dsDNA"]][,1],EnergyMap.AlphaCarbon[["Maps"]][["mT1.dsDNA"]][,2]),nlevels = 10,col = 'red',drawlabels = F,add=TRUE)
legend('topright',legend=c('mT1.Apo','mT1.dsDNA'),fill = c('black','red'),col = c('black','red'))

# save and cleanup

save(EnergyMap.AlphaCarbon,file = 'EnergyMap-AlphaCarbon.RData')
rm(list=ls())

#############################################################################

##### Ramachandran Clustering

####################################

# load data

rm(list=ls())

load('mT1-Apo_Input.RData')
mT1.Apo[['xyz']]=NULL

load('mT1-dsDNA_Input.RData')
mT1.dsDNA[['xyz']]=NULL

# transform

Ramachandran.Cluster.input.phi=rbind(rbind(mT1.Apo[['phi']][,,1],mT1.Apo[['phi']][,,2],mT1.Apo[['phi']][,,3],mT1.Apo[['phi']][,,4]),rbind(mT1.dsDNA[['phi']][,,1],mT1.dsDNA[['phi']][,,2],mT1.dsDNA[['phi']][,,3],mT1.dsDNA[['phi']][,,4]))[seq(1,2*prod(dim(mT1.Apo[['phi']])[c(1,3)]),2),231:460]
Ramachandran.Cluster.input.psi=rbind(rbind(mT1.Apo[['psi']][,,1],mT1.Apo[['psi']][,,2],mT1.Apo[['psi']][,,3],mT1.Apo[['psi']][,,4]),rbind(mT1.dsDNA[['psi']][,,1],mT1.dsDNA[['psi']][,,2],mT1.dsDNA[['psi']][,,3],mT1.dsDNA[['psi']][,,4]))[seq(1,2*prod(dim(mT1.Apo[['phi']])[c(1,3)]),2),231:460]
Ramachandran.Cluster.input.key=cbind(rep(c('mT1.Apo','mT1.dsDNA'),each=5001*4),rep(rep(1:4,each=5001),times=2))[seq(1,2*prod(dim(mT1.Apo[['phi']])[c(1,3)]),2),]
colnames(Ramachandran.Cluster.input.key)=c('SIM','RUN')
Ramachandran.Cluster.input.residues=paste(rep(c('B'),each=231),'_',mT1.Apo[['labels']][((mT1.Apo[['labels']][,5]=='protB') & mT1.Apo[['labels']][,2]=='CA'),4],'-',rep(4:234,times=1),sep='')[-c(231)]

Ramachandran.Cluster.input.x=cos(Ramachandran.Cluster.input.psi*pi/180)*sin(Ramachandran.Cluster.input.phi*pi/180)
Ramachandran.Cluster.input.y=sin(Ramachandran.Cluster.input.psi*pi/180)*sin(Ramachandran.Cluster.input.phi*pi/180)
Ramachandran.Cluster.input.z=cos(Ramachandran.Cluster.input.phi*pi/180)

Ramachandran.Cluster.input=list('x'=Ramachandran.Cluster.input.x,'y'=Ramachandran.Cluster.input.y,'z'=Ramachandran.Cluster.input.z,'phi'=Ramachandran.Cluster.input.phi,'psi'=Ramachandran.Cluster.input.psi,'key'=Ramachandran.Cluster.input.key,'residues'=Ramachandran.Cluster.input.residues)
rm(list=setdiff(ls(),'Ramachandran.Cluster.input'))

# analyze

library(dbscan)

Ramachandran.Cluster.dbscan.groups=matrix(-1,nrow=nrow(Ramachandran.Cluster.input[['phi']]),ncol = ncol(Ramachandran.Cluster.input[['phi']]))
for (i in 1:ncol(Ramachandran.Cluster.input[['phi']])){
  Ramachandran.Cluster.dbscan.groups[,i]=hdbscan(cbind(Ramachandran.Cluster.input[['x']][,i],Ramachandran.Cluster.input[['y']][,i],Ramachandran.Cluster.input[['z']][,i]),minPts = 200)[['cluster']]
}
Ramachandran.Cluster.dbscan=list('groups'=Ramachandran.Cluster.dbscan.groups,'key'=Ramachandran.Cluster.input[["key"]],'residues'=Ramachandran.Cluster.input[["residues"]])
rm(i,Ramachandran.Cluster.dbscan.groups)

# stats

ptest=function(SIM.1,SIM.2,RES){
  frame=matrix(0,ncol=2,nrow=(1+max(Ramachandran.Cluster.dbscan[['groups']][,RES])))
  for (j in 1:nrow(frame)){
    frame[j,1]=sum(Ramachandran.Cluster.dbscan[['groups']][Ramachandran.Cluster.dbscan[['key']][,1]==SIM.1,RES]==(j-1))
    frame[j,2]=sum(Ramachandran.Cluster.dbscan[['groups']][Ramachandran.Cluster.dbscan[['key']][,1]==SIM.2,RES]==(j-1))
  }
  p=chisq.test(frame[rowSums(frame)>0,])[['p.value']]
  return(p)
}

p.APOvDSDNA=rep(2,times=length(Ramachandran.Cluster.dbscan[['residues']]))
for(i in 1:length(p.APOvDSDNA)){
  p.APOvDSDNA[i]=ptest('mT1.Apo','mT1.dsDNA',i)
}
q.APOvDSDNA=p.adjust(p.APOvDSDNA[order(p.APOvDSDNA)],method = 'BH')
Stats.APOvDSDNA=cbind(Ramachandran.Cluster.dbscan[['residues']][order(p.APOvDSDNA)],p.APOvDSDNA[order(p.APOvDSDNA)],q.APOvDSDNA)
colnames(Stats.APOvDSDNA)=c('residues','p.value','q.value')

Ramachandran.Cluster.Stats=list('p.test'=ptest,'APOvDSDNA'=Stats.APOvDSDNA)

# reconcile and cleanup

Ramachandran.Cluster=list('input'=Ramachandran.Cluster.input,'dbscan'=Ramachandran.Cluster.dbscan,'stats'=Ramachandran.Cluster.Stats)

rm(list=setdiff(ls(),'Ramachandran.Cluster'))

# plot check

plot(1:230,as.numeric(Ramachandran.Cluster[['stats']][['APOvDSDNA']][match(Ramachandran.Cluster[["input"]][['residues']],Ramachandran.Cluster[['stats']][['APOvDSDNA']][,1]),3]==0),type = 'h',ylab='q < 0.05',xlab = 'Residues (N → C)',yaxt='n',xaxt='n',main = 'Ramachandran: mT1.Apo vs mT1.dsDNA')
axis(2,at = c(0,1),labels = c('no','yes'))
axis(1,at = 115.5,labels = 'Protomer-B')

# save and wipe

save(Ramachandran.Cluster,file='Ramachandran-Cluster.RData')
rm(list=ls())

#############################################################################

##### Cα Correlations

####################################

### mT1_Apo

# load and trim data

rm(list=ls())

load('mT1-Apo_Input.RData')
mT1.Apo.labels=mT1.Apo[['labels']]
mT1.Apo.xyz=mT1.Apo[['xyz']]
mT1.Apo.xyz.alphas=mT1.Apo.xyz[,((mT1.Apo.labels[,5]=='protA' | mT1.Apo.labels[,5]=='protB') & mT1.Apo.labels[,2]=='CA'),,]
residues=paste(rep(c('A','B'),each=231),'_',mT1.Apo[['labels']][((mT1.Apo[['labels']][,5]=='protA' | mT1.Apo[['labels']][,5]=='protB') & mT1.Apo[['labels']][,2]=='CA'),4],'-',rep(4:234,times=2),sep='')
rm(mT1.Apo,mT1.Apo.xyz,mT1.Apo.labels)

# reformat

Correlations.mT1.Apo.input=list('data'=array(mT1.Apo.xyz.alphas,dim=c(dim(mT1.Apo.xyz.alphas)[1:2],prod(dim(mT1.Apo.xyz.alphas)[3:4]))),'residues'=residues)
rm(list=setdiff(ls(),'Correlations.mT1.Apo.input'))

# analyze

correlate=function(A){
  temp.1=array(2,dim=c(dim(A)[2],dim(A)[2]))
  temp.2=array(apply(A,MARGIN = c(1,2),FUN = mean),dim=dim(A))
  temp.3=A-temp.2
  for (i in 1:dim(A)[2]){
    for (j in 1:dim(A)[2]){
      temp.1[i,j]=mean(colSums(temp.3[,i,]*temp.3[,j,]))/sqrt(mean(colSums(temp.3[,i,]^2))*mean(colSums(temp.3[,j,]^2)))
    }
  }
  return(temp.1)
}

Correlations.mT1.Apo.Cij=correlate(Correlations.mT1.Apo.input[['data']])

Correlations.mT1.Apo=list('input'=Correlations.mT1.Apo.input,'Cij'=Correlations.mT1.Apo.Cij)
Correlations.AlphaCarbon=list('mT1.Apo'=Correlations.mT1.Apo)
rm(list=setdiff(ls(),c('Correlations.AlphaCarbon')))

####################################

### mT1_4-mer

# load and trim data

load('mT1-4mer_Input.RData')
mT1.4mer.labels=mT1.4mer[['labels']]
mT1.4mer.xyz=mT1.4mer[['xyz']]
mT1.4mer.xyz.alphas=mT1.4mer.xyz[,((mT1.4mer.labels[,5]=='protA' | mT1.4mer.labels[,5]=='protB') & mT1.4mer.labels[,2]=='CA'),,]
residues=paste(rep(c('A','B'),each=231),'_',mT1.4mer[['labels']][((mT1.4mer[['labels']][,5]=='protA' | mT1.4mer[['labels']][,5]=='protB') & mT1.4mer[['labels']][,2]=='CA'),4],'-',rep(4:234,times=2),sep='')
rm(mT1.4mer,mT1.4mer.xyz,mT1.4mer.labels)

# reformat

Correlations.mT1.4mer.input=list('data'=array(mT1.4mer.xyz.alphas,dim=c(dim(mT1.4mer.xyz.alphas)[1:2],prod(dim(mT1.4mer.xyz.alphas)[3:4]))),'residues'=residues)
rm(list=setdiff(ls(),c('Correlations.AlphaCarbon','Correlations.mT1.4mer.input')))

# analyze

correlate=function(A){
  temp.1=array(2,dim=c(dim(A)[2],dim(A)[2]))
  temp.2=array(apply(A,MARGIN = c(1,2),FUN = mean),dim=dim(A))
  temp.3=A-temp.2
  for (i in 1:dim(A)[2]){
    for (j in 1:dim(A)[2]){
      temp.1[i,j]=mean(colSums(temp.3[,i,]*temp.3[,j,]))/sqrt(mean(colSums(temp.3[,i,]^2))*mean(colSums(temp.3[,j,]^2)))
    }
  }
  return(temp.1)
}

Correlations.mT1.4mer.Cij=correlate(Correlations.mT1.4mer.input[['data']])

Correlations.mT1.4mer=list('input'=Correlations.mT1.4mer.input,'Cij'=Correlations.mT1.4mer.Cij)
Correlations.AlphaCarbon[['mT1.4mer']]=Correlations.mT1.4mer
rm(list=setdiff(ls(),c('Correlations.AlphaCarbon')))

####################################

### FDR Stats

tests=array(2,dim = c(462,462,10))

for (i in 1:dim(tests)[3]){
  tests[,,i]=replicate(462^2,mean(runif(20001,-1,1)))
}

check.ref=seq(0,1,0.01)
props=array(0,dim = c(dim(tests)[3],length(check.ref)))
for (i in 1:dim(tests)[3]){
  for (j in 1:length(check.ref)){
    props[i,j]=sum(abs(tests[,,i])>=check.ref[j])/length(tests[,,i])
  }
}

Correlations.Stats=list('tests'=tests,'props'=props,'ref'=check.ref)
Correlations.AlphaCarbon[['Stats']]=Correlations.Stats

rm(list=setdiff(ls(),'Correlations.AlphaCarbon'))

# results check

heatmap(Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]],symm = TRUE,Colv=NA,Rowv = NA,RowSideColors = hcl.colors(sqrt(length(Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]])),palette = 'Red-Green'),col = hcl.colors(length(Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]]),palette = 'Red-Green'),zlim=c(-1,1),xlab = 'Residues (N → C)',ylab = 'Residues (N → C)',main='Cα Correlations: mT1-WT',labRow = NA,labCol = NA)

heatmap(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]],symm = TRUE,Colv=NA,Rowv = NA,RowSideColors = hcl.colors(sqrt(length(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]])),palette = 'Red-Green'),col = hcl.colors(length(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]]),palette = 'Red-Green'),zlim=c(-1,1),xlab = 'Residues (N → C)',ylab = 'Residues (N → C)',main='Cα Correlations: mT1-WT + 4-mer',labRow = NA,labCol = NA)

# save and cleanup

save(Correlations.AlphaCarbon,file='Correlations-AlphaCarbon.RData')
rm(list=ls())

#############################################################################

##### Figure Preparation

####################################

### Environment

rm(list=ls())
setwd('/media/wayne/DATA/MD_Simulations/ANALYSES/DNA-Binding_Project/')
setwd('~/Documents/Research/MD_Simulations/ANALYSES/DNA-Binding_Project/')

####################################

### Supplemental 1: Simulation RMSDs

# Load

rm(list=ls())
load('RMSD-Global.RData')

# S.1.1

plot(seq(0,1000,0.2),RMSD.Global[["mT1.Apo"]][["all"]][,1],type='l',lty='dotted',col='grey',main='mT1-WT',xlab = 'Time (ns)',ylab = 'Enzyme RMSD (Å)',ylim = c(0,max(RMSD.Global[["mT1.Apo"]][["all"]])))
lines(seq(0,1000,0.2),RMSD.Global[["mT1.Apo"]][["all"]][,2],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.Apo"]][["all"]][,3],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.Apo"]][["all"]][,4],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.Apo"]][["mean"]],type='l',lty='solid',col='black')

# S.1.1

plot(seq(0,1000,0.2),RMSD.Global[["mT1.4mer"]][["all"]][,1],type='l',lty='dotted',col='grey',main='mT1-WT + 4-mer',xlab = 'Time (ns)',ylab = 'Enzyme RMSD (Å)',ylim = c(0,max(RMSD.Global[["mT1.4mer"]][["all"]])))
lines(seq(0,1000,0.2),RMSD.Global[["mT1.4mer"]][["all"]][,2],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.4mer"]][["all"]][,3],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.4mer"]][["all"]][,4],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.4mer"]][["mean"]],type='l',lty='solid',col='black')

# S.1.1

plot(seq(0,1000,0.2),RMSD.Global[["mT1.dsDNA"]][["all"]][,1],type='l',lty='dotted',col='grey',main='mT1-WT + dsDNA',xlab = 'Time (ns)',ylab = 'Enzyme RMSD (Å)',ylim = c(0,max(RMSD.Global[["mT1.dsDNA"]][["all"]])))
lines(seq(0,1000,0.2),RMSD.Global[["mT1.dsDNA"]][["all"]][,2],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.dsDNA"]][["all"]][,3],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.dsDNA"]][["all"]][,4],type='l',lty='dotted',col='grey')
lines(seq(0,1000,0.2),RMSD.Global[["mT1.dsDNA"]][["mean"]],type='l',lty='solid',col='black')

####################################

### Figure 1: AlphaCarbon RMSFs

# Load

rm(list=ls())
load('RMSF-AlphaCarbon.RData')

# F.1.1

plot(1:231,RMSF.AlphaCarbon[["mT1.Apo"]][["mean"]],type='l',col='black',main='mTREX1 Backbone Fluctuations',ylab = 'RMSF (Å)',xlab='Residues (N → C)',ylim=c(0,50),xaxt='n')
lines(1:231,RMSF.AlphaCarbon[["mT1.dsDNA"]][["mean"]],col='red')
abline(v=c(160.5,172.5),col='green')
abline(v=c(19.5,24.5),col='purple')
abline(v=c(70.5,95.5),col='orange')
legend('topright',legend = c('mT1','mT1 + dsDNA'),fill = c('black','red'),col = c('black','red'))
axis(1,at = 115.5,labels = 'Protomer-B')

####################################

### Supplemental 2: AlphaCarbon Energy Map

# Load

rm(list=ls())
load('EnergyMap-AlphaCarbon.RData')
library(MASS)

# S.2.1

contour(kde2d(EnergyMap.AlphaCarbon[["Maps"]][["mT1.Apo"]][,1],EnergyMap.AlphaCarbon[["Maps"]][["mT1.Apo"]][,2]),nlevels = 10,col = 'black',drawlabels = F,xlim = c(-60,60),ylim = c(-50,50),main='Cα Dynamic Range (Protomer-B)',xlab='PC.1 (25.2%)',ylab='PC.2 (12.4%)')
contour(kde2d(EnergyMap.AlphaCarbon[["Maps"]][["mT1.dsDNA"]][,1],EnergyMap.AlphaCarbon[["Maps"]][["mT1.dsDNA"]][,2]),nlevels = 10,col = 'red',drawlabels = F,add=TRUE)
legend('topright',legend=c('mT1-WT','mT1-WT + dsDNA'),fill = c('black','red'),col = c('black','red'))

####################################

### Figure 2: Ramachandran Clustering

# Load

rm(list=ls())
load('Ramachandran-Cluster.RData')

# F.2.0

Ramachandran.Cluster[["input"]][['residues']][Ramachandran.Cluster[['stats']][['APOvDSDNA']][match(Ramachandran.Cluster[["input"]][['residues']],Ramachandran.Cluster[['stats']][['APOvDSDNA']][,1]),3]<=0.05]
((1:length(Ramachandran.Cluster[["input"]][['residues']]))+231)[Ramachandran.Cluster[['stats']][['APOvDSDNA']][match(Ramachandran.Cluster[["input"]][['residues']],Ramachandran.Cluster[['stats']][['APOvDSDNA']][,1]),3]<=0.05]


# F.2.1

plot(1:230,as.numeric(Ramachandran.Cluster[['stats']][['APOvDSDNA']][match(Ramachandran.Cluster[["input"]][['residues']],Ramachandran.Cluster[['stats']][['APOvDSDNA']][,1]),3]<=0.05),type = 'h',ylab='q < 0.05',xlab = 'Residues (N → C)',yaxt='n',xaxt='n',main = 'mT1-WT vs mT1-WT + dsDNA')
axis(2,at = c(0,1),labels = c('no','yes'))
axis(1,at = 115.5,labels = 'Protomer-B')

####################################

### Figure 3: Cα Correlations

# Load

rm(list=ls())
load('Correlations-AlphaCarbon.RData')

# F.3.1

heatmap(Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]],symm = TRUE,Colv=NA,Rowv = NA,RowSideColors = hcl.colors(sqrt(length(Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]])),palette = 'Red-Green'),col = hcl.colors(length(Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]]),palette = 'Red-Green'),zlim=c(-1,1),xlab = 'Residues (N → C)',ylab = 'Residues (N → C)',main='Cα Correlations: mT1-WT',labRow = NA,labCol = NA)

# F.3.2

heatmap(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]],symm = TRUE,Colv=NA,Rowv = NA,RowSideColors = hcl.colors(sqrt(length(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]])),palette = 'Red-Green'),col = hcl.colors(length(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]]),palette = 'Red-Green'),zlim=c(-1,1),xlab = 'Residues (N → C)',ylab = 'Residues (N → C)',main='Cα Correlations: mT1-WT + 4-mer',labRow = NA,labCol = NA)

# F.3.3

temp=(Correlations.AlphaCarbon[["mT1.4mer"]][["Cij"]]-Correlations.AlphaCarbon[["mT1.Apo"]][["Cij"]])
heatmap(temp,symm = TRUE,Colv=NA,Rowv = NA,RowSideColors = hcl.colors(sqrt(length(temp)),palette = 'Red-Green'),col = hcl.colors(length(temp),palette = 'Red-Green'),zlim=c(-1,1),xlab = 'Residues (N → C)',ylab = 'Residues (N → C)',main='Cα Correlations: Differentials',labRow = NA,labCol = NA)






























