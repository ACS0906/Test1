library(rrcov)
library(MASS)
library(robustbase)
library(robust)
library(pcaPP)
library(mvtnorm)
library(clusterGeneration)

#install.packages("mvtnorm")
#install.packages("pcaPP")
#install.packages("robustbase")
#install.packages("robustbase")
#install.packages("fit.models")

library(MASS)
D<-read.csv("4.3.csv",header=T)
#conta=matrix(rep(round(5*as.numeric(apply(D,2,max))+rnorm(9,0,1)),5),5)
DD=rbind(as.matrix(D),conta)
Y<-DD[,c(1,5,6)]
X<-DD[,c(2,4,8)]
nDD=cbind(X,Y)
dim(D)
v.mve<-cov2cor(CovMve(nDD)@cov) 
#v.bve<-cov2cor(MeanVar(nDD,Beta=0.1)$V)
r11<-v.mve[1:3,1:3]
r22<-v.mve[4:6,4:6]
r12<-v.mve[1:3,4:6]
r21<-t(r12)

###=====================================================   
#v.mcd<-CovMcd(nDD)@cov #?CovMcd
#v.ogk<-covOGK(nDD, sigmamu = scaleTau2)$wcov
###======================================================================
########### squre root matrix for r11
ev1<-eigen(r11) 
sev1<-sqrt(ev1$values)
dm1<-diag(1/sev1)
sr11<-ev1$vectors%*%dm1%*%t(ev1$vectors)
################# squre root matrix for r22
ev2<-eigen(r22) 
sev2<-sqrt(ev2$values)
dm2<-diag(1/sev2)
sr22<-ev2$vectors%*%dm2%*%t(ev2$vectors)
##################### Define Cross correlation matrix m
m<-sr11%*%r12%*%sr22
############## compute SVD OF m
sd<-svd(m)
##################### canonical correlation
can.corr<-sd$d
l=can.corr^2
################P.Value Calculation of canonical variate paires##########
Y<-DD[,c(1,5,6)]
X<-DD[,c(2,4,8)]
p=ncol(X);q=ncol(Y);n=nrow(X)
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
clx=matrix(rep(0,ncol(X)*ncol(Y)),ncol(X))
for(i in 1:ncol(Y))
{
L<-can.corr[i]*lx[,i]
clx[,i]<-L
}
clx<-clx
####################### cross loading for y

cly=matrix(rep(0,ncol(X)*ncol(Y)),ncol(Y))
for (i in 1:ncol(X))
{
L<-can.corr[i]*ly[,i]
cly[,i]<-L
}
#####################cross loading
clx<-clx
cly<-cly

################ propotion of Explained Variance for x set
pevx=rep(0,min(length(Y),length(X)))
for ( i in 1:min(length(Y),length(X))){
pex<-sum(lx[,i]^2)/length(lx[,i])
pevx[i]<-pex
}
### percent of variance for x set 

###########################propotion of Explained Variance for y set
pevy=rep(0,min(length(Y),length(X)))

for ( i in 1:min(length(Y),length(X))){
pey<-sum(ly[,i]^2)/length(ly[,i])
pevy[i]<-pey
}

### percent of variance  for y 
ppevx<-pevx*100
ppevy<-pevy*100

################Redundency for x
Rx=rep(0,min(p=3,q=3))
for ( i in 1:min(p=3,q=3)){
R1<-can.corr[i]^2*sum(lx[,i]^2)/length(lx[,i])
Rx[i]<-R1
}

#########################Redundency for y 
Ry=rep(0,min(p=3,q=3))
for ( i in 1:min(p=3,q=3)){
R2<-can.corr[i]^2*sum(ly[,i]^2)/length(ly[,i])
Ry[i]<-R2
}
Rx
Ry
Table=round(data.frame(Rx,Ry),3)
P.cca=cancor(X,Y)
##================================================================

#=====================================================
#Examples

data(hbk)

hbk.x <- data.matrix(hbk[, 1:3])

v.mcd<-attr(CovMcd(hbk.x),"cov")

## the following three statements are equivalent
c1 <- CovMcd(hbk.x, alpha = 0.75)
c2 <- CovMcd(hbk.x, control = CovControlMcd(alpha = 0.75))
## direct specification overrides control one:
c3 <- CovMcd(hbk.x, alpha = 0.75,
             control = CovControlMcd(alpha=0.95))
c1

###========================================================================