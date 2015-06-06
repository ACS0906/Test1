#======================Required Libraries==============================
library(MASS)
library(rrcov)
library(MASS)
library(robustbase)
library(robust)
library(pcaPP)
library(mvtnorm)
library(clusterGeneration)
#=====================Data Generation =============================
set.seed(0001)
mu=rep(0,6)
v=genPositiveDefMat("unifcorrmat",dim=6)$Sigma
Data=mvrnorm(n=100,mu,v);colnames(Data)=c("x1","x2","x3","y1","y2","y3")
#========================Correlation Matrix=========================
Orgi=cor(Data)
r11<-Orgi[1:6,1:6]
r22<-Orgi[7:8,7:8]
r12<-Orgi[1:6,7:8]
r21<-t(r12)

#===================== squre root matrix for r11===================
ev1<-eigen(r11) 
sev1<-sqrt(ev1$values)
dm1<-diag(1/sev1)
sr11<-ev1$vectors%*%dm1%*%t(ev1$vectors)
#===================== squre root matrix for r22====================
ev2<-eigen(r22) 
sev2<-sqrt(ev2$values)
dm2<-diag(1/sev2)
sr22<-ev2$vectors%*%dm2%*%t(ev2$vectors)
#====================== Define Cross correlation matrix m ==========
m<-sr11%*%r12%*%sr22
#====================== compute SVD OF m ===========================
sd<-svd(m)
#====================== canonical correlation========================
can.corr<-sd$d   #Canonical Correlations
l=can.corr^2
################P.Value Calculation of canonical variate paires##########
p=ncol(r11);q=ncol(r22);n=dim(Data)[1]
L=PV=df=matrix(rep(0,min(p,q)),nrow=min(p,q))
 
for(m in 1:min(p,q))
{
 df[m]<-(p-m+1)*(q-m+1)
 df  
 X=0
 for(k in m:min(p,q))
    { 
     X=X+log(1-l[k])
     }
    L[m]=-((n-1)-.5*(p+q+1))*X
    PV=pchisq(L,df=df,lower.tail=F) 
}
orginal=L;PV;df
Table=round(data.frame(CanCorr=can.corr,Chisquare=L,df=df,p.value=PV),3)

################## canonical weights 
a<-sr11%*%sd$u
a[,colSums(a)<0]=-a[,colSums(a)<0]
b<-sr22%*%sd$v
b[,colSums(b)<0]=-b[,colSums(b)<0]
################### caninical loading 
lx<-r11%*%a
ly<-r22%*%b
#hist(lx[,1])
######################## cross loading for x
clx=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r11))
for(i in 1:ncol(r22))
{
L<-can.corr[i]*lx[,i]
clx[,i]<-L
}
clx<-clx
####################### cross loading for y
cly=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r22))
for (i in 1:ncol(r11))
{
L<-can.corr[i]*ly[,i]
cly[,i]<-L
}
#####################cross loading
clx<-clx
cly<-cly

################ propotion of Explained Variance for x set
pevx=rep(0,min(ncol(r22),ncol(r11)))
for ( i in 1:min(ncol(r22),ncol(r11))){
pex<-sum(lx[,i]^2)/length(lx[,i])
pevx[i]<-pex
}
### percent of variance for x set 
###########################propotion of Explained Variance for y set
pevy=rep(0,min(ncol(r22),ncol(r11)))

for ( i in 1:min(ncol(r22),ncol(r11))){
pey<-sum(ly[,i]^2)/length(ly[,i])
pevy[i]<-pey
}
### percent of variance  for y 
ppevx<-pevx*100
ppevy<-pevy*100

################Redundency for x
Rx=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R1<-can.corr[i]^2*sum(lx[,i]^2)/length(lx[,i])
Rx[i]<-R1
}

#########################Redundency for y 
Ry=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R2<-can.corr[i]^2*sum(ly[,i]^2)/length(ly[,i])
Ry[i]<-R2
}
Table=round(data.frame(Rx,Ry),3)


#===================================================================
##==============Different Covariance matrix========================
#===================================================================
v.bve<-cov2cor(MeanVar(Data,Beta=0.1)$V)
v.mve<-cov2cor(CovMve(Data)@cov) 
v.mcd<-cov2cor(CovMcd(Data)@cov) #?CovMcd
v.ogk<-cov2cor(covOGK(Data, sigmamu = scaleTau2)$wcov)


#=======================================================================
##=============Robustification by Beta divergence ================
##======================================================================
r11<-v.bve[1:6,1:6]
r22<-v.bve[7:8,7:8]
r12<-v.bve[1:6,7:8]
r21<-t(r12)
#===================== squre root matrix for r11===================
ev1<-eigen(r11) 
sev1<-sqrt(ev1$values)
dm1<-diag(1/sev1)
sr11<-ev1$vectors%*%dm1%*%t(ev1$vectors)
#===================== squre root matrix for r22====================
ev2<-eigen(r22) 
sev2<-sqrt(ev2$values)
dm2<-diag(1/sev2)
sr22<-ev2$vectors%*%dm2%*%t(ev2$vectors)
#====================== Define Cross correlation matrix m ==========
m<-sr11%*%r12%*%sr22
#====================== compute SVD OF m ===========================
sd<-svd(m)
#====================== canonical correlation========================
can.corr<-sd$d   #Canonical Correlations
l=can.corr^2
################P.Value Calculation of canonical variate paires##########
p=ncol(r11);q=ncol(r22);n=dim(Data)[1]
L=PV=df=matrix(rep(0,min(p,q)),nrow=min(p,q))
 
for(m in 1:min(p,q))
{
 df[m]<-(p-m+1)*(q-m+1)
 df  
 X=0
 for(k in m:min(p,q))
    { 
     X=X+log(1-l[k])
     }
    L[m]=-((n-1)-.5*(p+q+1))*X
    PV=pchisq(L,df=df,lower.tail=F) 
}
L;PV;df
Table1=round(data.frame(CanCorr=can.corr,Chisquare=L,df=df,p.value=PV),3)

################## canonical weights 
a<-sr11%*%sd$u
a[,colSums(a)<0]=-a[,colSums(a)<0]
b<-sr22%*%sd$v
b[,colSums(b)<0]=-b[,colSums(b)<0]
################### caninical loading 
lx<-r11%*%a
ly<-r22%*%b
#hist(lx[,1])
######################## cross loading for x
clx=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r11))
for(i in 1:ncol(r22))
{
L<-can.corr[i]*lx[,i]
clx[,i]<-L
}
clx<-clx
####################### cross loading for y
cly=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r22))
for (i in 1:ncol(r11))
{
L<-can.corr[i]*ly[,i]
cly[,i]<-L
}
#####################cross loading
clx<-clx
cly<-cly

################ propotion of Explained Variance for x set
pevx=rep(0,min(ncol(r22),ncol(r11)))
for ( i in 1:min(ncol(r22),ncol(r11))){
pex<-sum(lx[,i]^2)/length(lx[,i])
pevx[i]<-pex
}
### percent of variance for x set 
###########################propotion of Explained Variance for y set
pevy=rep(0,min(ncol(r22),ncol(r11)))

for ( i in 1:min(ncol(r22),ncol(r11))){
pey<-sum(ly[,i]^2)/length(ly[,i])
pevy[i]<-pey
}
### percent of variance  for y 
ppevx<-pevx*100
ppevy<-pevy*100

################Redundency for x
Rx=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R1<-can.corr[i]^2*sum(lx[,i]^2)/length(lx[,i])
Rx[i]<-R1
}

#########################Redundency for y 
Ry=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R2<-can.corr[i]^2*sum(ly[,i]^2)/length(ly[,i])
Ry[i]<-R2
}
Table=round(data.frame(Rx,Ry),3)




#=============================================================================
##=============Robustification by Mve estimator ================
##===========================================================================

r11<-v.mve[1:3,1:3]
r22<-v.mve[4:6,4:6]
r12<-v.mve[1:3,4:6]
r21<-t(r12)
#===================== squre root matrix for r11===================
ev1<-eigen(r11) 
sev1<-sqrt(ev1$values)
dm1<-diag(1/sev1)
sr11<-ev1$vectors%*%dm1%*%t(ev1$vectors)
#===================== squre root matrix for r22====================
ev2<-eigen(r22) 
sev2<-sqrt(ev2$values)
dm2<-diag(1/sev2)
sr22<-ev2$vectors%*%dm2%*%t(ev2$vectors)
#====================== Define Cross correlation matrix m ==========
m<-sr11%*%r12%*%sr22
#====================== compute SVD OF m ===========================
sd<-svd(m)
#====================== canonical correlation========================
can.corr<-sd$d   #Canonical Correlations
l=can.corr^2
################P.Value Calculation of canonical variate paires##########
p=ncol(r11);q=ncol(r22);n=dim(Data)[1]
L=PV=df=matrix(rep(0,min(p,q)),nrow=min(p,q))
 
for(m in 1:min(p,q))
{
 df[m]<-(p-m+1)*(q-m+1)
 df  
 X=0
 for(k in m:min(p,q))
    { 
     X=X+log(1-l[k])
     }
    L[m]=-((n-1)-.5*(p+q+1))*X
    PV=pchisq(L,df=df,lower.tail=F) 
}
L;PV;df
Table=round(data.frame(CanCorr=can.corr,Chisquare=L,df=df,p.value=PV),3)

################## canonical weights 
a<-sr11%*%sd$u
a[,colSums(a)<0]=-a[,colSums(a)<0]
b<-sr22%*%sd$v
b[,colSums(b)<0]=-b[,colSums(b)<0]
################### caninical loading 
lx<-r11%*%a
ly<-r22%*%b
#hist(lx[,1])
######################## cross loading for x
clx=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r11))
for(i in 1:ncol(r22))
{
L<-can.corr[i]*lx[,i]
clx[,i]<-L
}
clx<-clx
####################### cross loading for y
cly=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r22))
for (i in 1:ncol(r11))
{
L<-can.corr[i]*ly[,i]
cly[,i]<-L
}
#####################cross loading
clx<-clx
cly<-cly

################ propotion of Explained Variance for x set
pevx=rep(0,min(ncol(r22),ncol(r11)))
for ( i in 1:min(ncol(r22),ncol(r11))){
pex<-sum(lx[,i]^2)/length(lx[,i])
pevx[i]<-pex
}
### percent of variance for x set 
###########################propotion of Explained Variance for y set
pevy=rep(0,min(ncol(r22),ncol(r11)))

for ( i in 1:min(ncol(r22),ncol(r11))){
pey<-sum(ly[,i]^2)/length(ly[,i])
pevy[i]<-pey
}
### percent of variance  for y 
ppevx<-pevx*100
ppevy<-pevy*100

################Redundency for x
Rx=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R1<-can.corr[i]^2*sum(lx[,i]^2)/length(lx[,i])
Rx[i]<-R1
}

#########################Redundency for y 
Ry=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R2<-can.corr[i]^2*sum(ly[,i]^2)/length(ly[,i])
Ry[i]<-R2
}
Table=round(data.frame(Rx,Ry),3)

#=============================================================================
##=============Robustification by MCD  estimator ================
##===========================================================================

r11<-v.mcd[1:3,1:3]
r22<-v.mcd[4:6,4:6]
r12<-v.mcd[1:3,4:6]
r21<-t(r12)
#===================== squre root matrix for r11===================
ev1<-eigen(r11) 
sev1<-sqrt(ev1$values)
dm1<-diag(1/sev1)
sr11<-ev1$vectors%*%dm1%*%t(ev1$vectors)
#===================== squre root matrix for r22====================
ev2<-eigen(r22) 
sev2<-sqrt(ev2$values)
dm2<-diag(1/sev2)
sr22<-ev2$vectors%*%dm2%*%t(ev2$vectors)
#====================== Define Cross correlation matrix m ==========
m<-sr11%*%r12%*%sr22
#====================== compute SVD OF m ===========================
sd<-svd(m)
#====================== canonical correlation========================
can.corr<-sd$d   #Canonical Correlations
l=can.corr^2
################P.Value Calculation of canonical variate paires##########
p=ncol(r11);q=ncol(r22);n=dim(Data)[1]
L=PV=df=matrix(rep(0,min(p,q)),nrow=min(p,q))
 
for(m in 1:min(p,q))
{
 df[m]<-(p-m+1)*(q-m+1)
 df  
 X=0
 for(k in m:min(p,q))
    { 
     X=X+log(1-l[k])
     }
    L[m]=-((n-1)-.5*(p+q+1))*X
    PV=pchisq(L,df=df,lower.tail=F) 
}
L;PV;df
Table=round(data.frame(CanCorr=can.corr,Chisquare=L,df=df,p.value=PV),3)

################## canonical weights 
a<-sr11%*%sd$u
a[,colSums(a)<0]=-a[,colSums(a)<0]
b<-sr22%*%sd$v
b[,colSums(b)<0]=-b[,colSums(b)<0]
################### caninical loading 
lx<-r11%*%a
ly<-r22%*%b
#hist(lx[,1])
######################## cross loading for x
clx=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r11))
for(i in 1:ncol(r22))
{
L<-can.corr[i]*lx[,i]
clx[,i]<-L
}
clx<-clx
####################### cross loading for y
cly=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r22))
for (i in 1:ncol(r11))
{
L<-can.corr[i]*ly[,i]
cly[,i]<-L
}
#####################cross loading
clx<-clx
cly<-cly

################ propotion of Explained Variance for x set
pevx=rep(0,min(ncol(r22),ncol(r11)))
for ( i in 1:min(ncol(r22),ncol(r11))){
pex<-sum(lx[,i]^2)/length(lx[,i])
pevx[i]<-pex
}
### percent of variance for x set 
###########################propotion of Explained Variance for y set
pevy=rep(0,min(ncol(r22),ncol(r11)))

for ( i in 1:min(ncol(r22),ncol(r11))){
pey<-sum(ly[,i]^2)/length(ly[,i])
pevy[i]<-pey
}
### percent of variance  for y 
ppevx<-pevx*100
ppevy<-pevy*100

################Redundency for x
Rx=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R1<-can.corr[i]^2*sum(lx[,i]^2)/length(lx[,i])
Rx[i]<-R1
}

#########################Redundency for y 
Ry=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R2<-can.corr[i]^2*sum(ly[,i]^2)/length(ly[,i])
Ry[i]<-R2
}
Table=round(data.frame(Rx,Ry),3)

#=========================================================================
#===========Robustification by OGK estimator=================
##========================================================================
r11<-v.ogk[1:3,1:3]
r22<-v.ogk[4:6,4:6]
r12<-v.ogk[1:3,4:6]
r21<-t(r12)
#===================== squre root matrix for r11===================
ev1<-eigen(r11) 
sev1<-sqrt(ev1$values)
dm1<-diag(1/sev1)
sr11<-ev1$vectors%*%dm1%*%t(ev1$vectors)
#===================== squre root matrix for r22====================
ev2<-eigen(r22) 
sev2<-sqrt(ev2$values)
dm2<-diag(1/sev2)
sr22<-ev2$vectors%*%dm2%*%t(ev2$vectors)
#====================== Define Cross correlation matrix m ==========
m<-sr11%*%r12%*%sr22
#====================== compute SVD OF m ===========================
sd<-svd(m)
#====================== canonical correlation========================
can.corr<-sd$d   #Canonical Correlations
l=can.corr^2
################P.Value Calculation of canonical variate paires##########
p=ncol(r11);q=ncol(r22);n=dim(Data)[1]
L=PV=df=matrix(rep(0,min(p,q)),nrow=min(p,q))
 
for(m in 1:min(p,q))
{
 df[m]<-(p-m+1)*(q-m+1)
 df  
 X=0
 for(k in m:min(p,q))
    { 
     X=X+log(1-l[k])
     }
    L[m]=-((n-1)-.5*(p+q+1))*X
    PV=pchisq(L,df=df,lower.tail=F) 
}
L;PV;df
Table=round(data.frame(CanCorr=can.corr,Chisquare=L,df=df,p.value=PV),3)

################## canonical weights 
a<-sr11%*%sd$u
a[,colSums(a)<0]=-a[,colSums(a)<0]
b<-sr22%*%sd$v
b[,colSums(b)<0]=-b[,colSums(b)<0]
################### caninical loading 
lx<-r11%*%a
ly<-r22%*%b
#hist(lx[,1])
######################## cross loading for x
clx=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r11))
for(i in 1:ncol(r22))
{
L<-can.corr[i]*lx[,i]
clx[,i]<-L
}
clx<-clx
####################### cross loading for y
cly=matrix(rep(0,ncol(r11)*ncol(r22)),ncol(r22))
for (i in 1:ncol(r11))
{
L<-can.corr[i]*ly[,i]
cly[,i]<-L
}
#####################cross loading
clx<-clx
cly<-cly

################ propotion of Explained Variance for x set
pevx=rep(0,min(ncol(r22),ncol(r11)))
for ( i in 1:min(ncol(r22),ncol(r11))){
pex<-sum(lx[,i]^2)/length(lx[,i])
pevx[i]<-pex
}
### percent of variance for x set 
###########################propotion of Explained Variance for y set
pevy=rep(0,min(ncol(r22),ncol(r11)))

for ( i in 1:min(ncol(r22),ncol(r11))){
pey<-sum(ly[,i]^2)/length(ly[,i])
pevy[i]<-pey
}
### percent of variance  for y 
ppevx<-pevx*100
ppevy<-pevy*100

################Redundency for x
Rx=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R1<-can.corr[i]^2*sum(lx[,i]^2)/length(lx[,i])
Rx[i]<-R1
}

#########################Redundency for y 
Ry=rep(0,min(p,q))
for ( i in 1:min(p,q)){
R2<-can.corr[i]^2*sum(ly[,i]^2)/length(ly[,i])
Ry[i]<-R2
}
Table=round(data.frame(Rx,Ry),3)

##  Classical CCA End  === We observe same result for all methods









