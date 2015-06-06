

set.seed(0001)
mu=rep(0,3)
v=genPositiveDefMat("unifcorrmat",dim=3)$Sigma
Data=mvrnorm(n=10,mu,v);colnames(Data)=c("x1","x2","x3")
conta=matrix(rep(round(3*as.numeric(apply(Data,2,max))+rnorm(3,0,1)),3),3)


nn=6
rOut=4;cOut=4;nv=3
#=============================Outlier Function===========
Outlier<-function(Data,rOut,cOut,nv)
{

OutIndex=sample(1:nn)[1:round(rOut*nn/100)] 
   DataOut=Data 
   for (i1 in 1:length(OutIndex))
   {
    OutCol=sample(1:nv)[1:round(cOut*nv/100)] #[1:(cOut*round(nv/4))]
    OUT=as.matrix(DataOut[,OutCol])
    OUT[OutIndex[i1],]<-2.5*OUT[OutIndex[i1],]
    DataOut[,OutCol]<-OUT  
    } 
    if(rOut==0){DataOut=Data}
return(DataOut)
}

Outlier(Data,rOut,cOut,nv)
# sample()'s surprise -- example
x <- 1:10
    sample(x[x >  8]) # length 2
    sample(x[x >  9]) # oops -- length 10!
    sample(x[x > 10]) # length 0

resample <- function(x, ...) x[sample.int(length(x), ...)]
resample(x[x >  8]) # length 2
resample(x[x >  9]) # length 1
resample(x[x > 10]) # length 0

