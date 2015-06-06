
set.seed(0001)
mu=rep(0,10)
v=genPositiveDefMat("unifcorrmat",dim=10)$Sigma
Data=mvrnorm(n=100,mu,v)#;colnames(Data)=c("x1","x2","x3","x4","x5","y1","y2","y3","y4","y5")

set.seed(01)
mu=rep(0,20)
v=genPositiveDefMat("unifcorrmat",dim=20)$Sigma
Data=mvrnorm(n=100,mu,v);colnames(Data)=c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","y1","y2","y3","y4","y5","y6","y7","y8","y9","y10")

set.seed(01)
mu=rep(0,50)
v=genPositiveDefMat("unifcorrmat",dim=50)$Sigma
Data=mvrnorm(n=90,mu,v);colnames(Data)=c(1:50)

set.seed(0001)
mu=rep(0,100)
v=genPositiveDefMat("unifcorrmat",dim=100)$Sigma
conta=mvrnorm(n=100,mu,v);colnames=c(1:100)

Datav=cbind(Data,conta)

dim(Datav)

Orgi=cor(Datav)
r11<-Orgi[1:3,1:3]
r22<-Orgi[4:6,4:6]
r12<-Orgi[1:3,4:6]
r21<-t(r12)



set.seed(0001)
mu=rep(0,6)
v=genPositiveDefMat("unifcorrmat",dim=6)$Sigma
Data=mvrnorm(n=100,mu,v)
conta=mvrnorm(100,rep(0,49),genPositiveDefMat("unifcorrmat",dim=49)$Sigma)
Data=cbind(Data,conta);dim(Data)



##==============Different Covariance matrix========================
#===================================================================
# We observe that there also same result though the data are contaminated
##========================================================================

v.bve<-cov2cor(MeanVar(Data,Beta=0.1)$V)
v.mve<-cov2cor(CovMve(Datav)@cov) 
v.mcd<-cov2cor(CovMcd(Datav)@cov) #?CovMcd
v.ogk<-cov2cor(covOGK(Data, sigmamu = scaleTau2)$wcov)

#=======================================================================
##=============Robustification by Beta divergence ================
##======================================================================
r11<-v.bve[1:3,1:3]
r22<-v.bve[4:6,4:6]
r12<-v.bve[1:3,4:6]
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
can.corr1<-sd$d;can.corr1
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
#L;PV;df
Table1=round(data.frame(CanCorr=can.corr1,Chisquare=L,df=df,p.value=PV),3);Table1


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
can.corr2<-sd$d;can.corr2
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
#L;PV;df
Table1=round(data.frame(CanCorr=can.corr2,Chisquare=L,df=df,p.value=PV),3);Table1

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
can.corr3<-sd$d;can.corr3
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
#L;PV;df
Table1=round(data.frame(CanCorr=can.corr3,Chisquare=L,df=df,p.value=PV),3);Table1



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
can.corr4<-sd$d ;can.corr4
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
#L;PV;df
Table1=round(data.frame(CanCorr=can.corr4,Chisquare=L,df=df,p.value=PV),3);Table1
###return(list(Bcca=cancorr1,MVE=cancorr2,MCD=cancorr3,OGK=cancorr4))

