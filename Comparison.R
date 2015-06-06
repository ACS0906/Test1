####====Requared Libraries ===========================
library(rrcov)
library(MASS)
library(robustbase)
library(robust)
library(pcaPP)
library(mvtnorm)
library(clusterGeneration)
#install.packages("clusterGeneration")
## =========Generating Data matrix ===============================
mux=muy=rep(0,3)
vx=genPositiveDefMat("unifcorrmat",dim=3)$Sigma #genPositiveDefMat("unifcorrmat",dim=3)
vy=genPositiveDefMat("eigen",dim=3)$Sigma
X=mvrnorm(n=100,mux,vx);colnames(X)=c("x1","x2","x3");cor(X)
Y=mvrnorm(n=100,mux,vy);colnames(Y)=c("y1","y2","y3");cor(Y)
Z=cbind(X,Y)
A=matrix(rnorm(36,6,.5),6,6);colnames(A)=c("x1","x2","x3","y1","y2","y3")
ZZ=Z%*%A
#pairs(ZZ)
#plot(ZZ)
#pairs(Z)
#===================================================================
##==============Different Covariance matrix========================
#===================================================================
v.bve<-cov2cor(MeanVar(ZZ,Beta=0.1)$V)
v.mve<-cov2cor(CovMve(ZZ)@cov) 
v.mcd<-cov2cor(CovMcd(ZZ)@cov) #?CovMcd
v.ogk<-cov2cor(covOGK(ZZ, sigmamu = scaleTau2)$wcov)

####==========CCA calculaton===========================================

r11<-v.mve[1:3,1:3]
r22<-v.mve[4:6,4:6]
r12<-v.mve[1:3,4:6]
r21<-t(r12)









Examples
genPositiveDefMat("unifcorrmat",dim=4)
Examples

Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
var(mvrnorm(n=1000, rep(0, 2), Sigma))
var(mvrnorm(n=1000, rep(0, 2), Sigma, empirical = TRUE))
