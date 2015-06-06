##### CCA Calculation with outlier ##############
#================================================
source(MeanVar.r)
library(MASS)
Data=rbind(matrix(rnorm(90,2,1),ncol=6),matrix(rnorm(30,10,1),ncol=6))
colnames(Data)=c("x1","x2","x3","x4","x5","x6")
plot(Data) #cor(Data)
x=Data[,c(1,2,3)]
y=Data[,c(4,5,6)]
Rxx=cor(x)
Ryy=cor(y)
Rxy=cor(x,y)
Ryx=cor(y,x)
CCA=function(Rxx,Ryy,Rxy,Ryx)
{
  Px<-solve(Rxx)%*%Rxy%*%solve(Ryy)%*%t(Rxy);Px
  Py<-solve(Ryy)%*%t(Rxy)%*%solve(Rxx)%*%Rxy;Py
  Egn1=eigen(Px);Egn2=eigen(Py)
  l=eigen(Px)$values
  #eigen(Px)
  a<-eigen(Px)$vectors#[,1]
  b<-eigen(Py)$vectors#[,1]
  ccc=sqrt(l)
  p<-length(Y);p
  q<-length(X);q
  n<-nrow(X);n
  CHI<- -((n-1)-.5*(p+q+1))*log(1-l)
  Rxxi=Rxx%*%a    ## Caninical Loading :Measure of simple correlation of original variable
  Ryyi=Ryy%*%b     ##with corresponding canonical variable.
  
  rxy=Ryyi%*%sqrt(l) #Cross loading
  ryx=Rxxi%*%sqrt(l)
  #Rsqiy=
  Rsqygx=(l/p)%*%t(Ryyi)%*%Ryyi # Proportion of explained variance
  Rsqxgy=(l/q)%*%t(Rxxi)%*%Rxx

  return(list(Canonical=ccc,Chisquare=CHI,LoadingX=Rxxi,
  LoadingY=Ryyi,CrossloadXY=rxy,CrosssloadYX=ryx))
}
Result=CCA(Rxx,Ryy,Rxy,Ryx)
RR=cancor(x,y)

###========================================================
##===ROBUST CCA CALCULATION WITH OUTLIER ===############################
DataR=cov2cor(MeanVar(Data,Beta=0.1)$V)
Rxxr=DataR[1:3,1:3]
Ryyr=DataR[4:6,4:6]
Rxyr=DataR[1:3,4:6]
Ryxr=DataR[4:6,1:3]
CCA=function(Rxxr,Ryyr,Rxyr,Ryxr)
{
  
  Pxr<-solve(Rxxr)%*%Rxyr%*%solve(Ryyr)%*%t(Rxyr);Pxr
  Pyr<-solve(Ryyr)%*%t(Rxyr)%*%solve(Rxxr)%*%Rxyr;Pyr
  Egn1=eigen(Pxr);Egn2=eigen(Pyr)
  lr=eigen(Pxr)$values
  #eigen(Pxr)
  ar<-eigen(Pxr)$vectors#[,1]
  br<-eigen(Pyr)$vectors#[,1]
  rccc=sqrt(l)
  p<-length(Y);p
  q<-length(X);q
  n<-nrow(Data);n
  CHIr<- -((n-1)-.5*(p+q+1))*log(1-l)
  Rxxir=Rxxr%*%ar    ## Caninical Loading :Measure of simple correlation of original variable
  Ryyir=Ryyr%*%br     ##with corresponding canonical variable.
  
  rxyr=Ryyir%*%sqrt(lr) #Cross loading
  ryxr=Rxxir%*%sqrt(lr)
  #Rsqiy=
  Rsqygx=(l/p)%*%t(Ryyir)%*%Ryyir # Proportion of explained variance
  Rsqxgy=(l/q)%*%t(Rxxir)%*%Rxxr

  return(list(Canonical=rccc,Chisquare=CHIr,LoadingXr=Rxxir,
  LoadingYr=Ryyir,CrossloadXYr=rxyr,CrosssloadYXr=ryxr))
}
Resultr=CCA(Rxxr,Ryyr,Rxyr,Ryxr)


###======================================================================
###Without outlier CCA========================================
DD=matrix(rnorm(120,2,1),ncol=6)
colnames(DD)=c("x1","x2","x3","x4","x5","x6")
plot(DD) #cor(Data)
x=DD[,c(1,2,3)]
y=DD[,c(4,5,6)]
Rxx=cor(x)
Ryy=cor(y)
Rxy=cor(x,y)
Ryx=cor(y,x)
CCA=function(Rxx,Ryy,Rxy,Ryx)
{
  Px<-solve(Rxx)%*%Rxy%*%solve(Ryy)%*%t(Rxy);Px
  Py<-solve(Ryy)%*%t(Rxy)%*%solve(Rxx)%*%Rxy;Py
  Egn1=eigen(Px);Egn2=eigen(Py)
  l=eigen(Px)$values
  #eigen(Px)
  a<-eigen(Px)$vectors#[,1]
  b<-eigen(Py)$vectors#[,1]
  ccc=sqrt(l)
  p<-length(Y);p
  q<-length(X);q
  n<-nrow(X);n
  CHI<- -((n-1)-.5*(p+q+1))*log(1-l)
  Rxxi=Rxx%*%a    ## Caninical Loading :Measure of simple correlation of original variable
  Ryyi=Ryy%*%b     ##with corresponding canonical variable.
  
  rxy=Ryyi%*%sqrt(l) #Cross loading
  ryx=Rxxi%*%sqrt(l)
  #Rsqiy=
  Rsqygx=(l/p)%*%t(Ryyi)%*%Ryyi # Proportion of explained variance
  Rsqxgy=(l/q)%*%t(Rxxi)%*%Rxx

  return(list(Canonical=ccc,Chisquare=CHI,LoadingX=Rxxi,
  LoadingY=Ryyi,CrossloadXY=rxy,CrosssloadYX=ryx))
}
Result=CCA(Rxx,Ryy,Rxy,Ryx)
Resultp=cancor(x,y)


###========================================================
##===ROBUST CCA CALCULATION ===############################
DataR=cov2cor(MeanVar(DD,Beta=0.1)$V)
Rxxr=DataR[1:3,1:3]
Ryyr=DataR[4:6,4:6]
Rxyr=DataR[1:3,4:6]
Ryxr=DataR[4:6,1:3]
CCA=function(Rxxr,Ryyr,Rxyr,Ryxr)
{ 
  Pxr<-solve(Rxxr)%*%Rxyr%*%solve(Ryyr)%*%t(Rxyr);Pxr
  Pyr<-solve(Ryyr)%*%t(Rxyr)%*%solve(Rxxr)%*%Rxyr;Pyr
  Egn1=eigen(Pxr);Egn2=eigen(Pyr)
  lr=eigen(Pxr)$values
  #eigen(Pxr)
  ar<-eigen(Pxr)$vectors#[,1]
  br<-eigen(Pyr)$vectors#[,1]
  rccc=sqrt(l)
  p<-length(Y);p
  q<-length(X);q
  n<-nrow(Data);n
  CHIr<- -((n-1)-.5*(p+q+1))*log(1-l)
  Rxxir=Rxxr%*%ar    ## Caninical Loading :Measure of simple correlation of original variable
  Ryyir=Ryyr%*%br     ##with corresponding canonical variable.
  
  rxyr=Ryyir%*%sqrt(lr) #Cross loading
  ryxr=Rxxir%*%sqrt(lr)
  #Rsqiy=
  Rsqygx=(l/p)%*%t(Ryyir)%*%Ryyir # Proportion of explained variance
  Rsqxgy=(l/q)%*%t(Rxxir)%*%Rxxr
  return(list(Canonical=rccc,Chisquare=CHIr,LoadingXr=Rxxir,
  LoadingYr=Ryyir,CrossloadXYr=rxyr,CrosssloadYXr=ryxr))
}
Resultr=CCA(Rxxr,Ryyr,Rxyr,Ryxr)


###======================================================================
