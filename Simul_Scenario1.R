
rm(list=ls(all=TRUE))

source("Moment_SMSNT.R")
source("Integral_nu_float.R")
source("algEM.fmr.ST.crFinal.R")

#install.packages("mixsmsn")
library(mixsmsn)

################################################################################
# Simulation Studies
# Simulation 1
################################################################################

n<-1000 # n : 100, 200, 500, 1000

x<-cbind(1, runif(n,1,5),runif(n,0,1),rbinom(n, 1, 0.6))

y<-matrix(0,n,1)
gru<-matrix(0,n,1)

betas1<-c(1,-1,-4,1)
Sigma1 <- 1
mu1 <- x%*%betas1
shape1<- -2

betas2<-c(1,2,3,1)
Sigma2 <- 2
mu2 <- x%*%betas2
shape2<- 3

nu<-4 # In the simulation study, nu is fixed

pii<-c(0.7,0.3) # mixing proportions

perc<-0.05 # censoring level = 0.05, 0.1, 0.2, 0.5

group.data<-list()

cutof<-matrix(0, nrow=100, ncol=length(perc))

################################################################################
## calculate a censored point corresponds to a censoring level
################################################################################

for(j in 1:100){ 
  for (i in 1:n){
    arg1 = list(mu=mu1[i], sigma2=Sigma1, shape=shape1, nu=nu)
    arg2 = list(mu=mu2[i], sigma2=Sigma2, shape=shape2, nu=nu)
    compl= rmix(1, pii, "Skew.t", list(arg1,arg2), cluster=TRUE) # change to ("t", "Skew.t", "Skew.cn", "Skew.slash",
    #"Skew.normal", "Normal" from the mixsmsns package
    y[i] = compl$y
    gru[i] = compl$cluster
  }
  
  ## creating censored observations (left censoring)
    aa=sort(y,decreasing=FALSE)
    cutof[j]<-aa[ceiling(perc*n)]
}
cutoff <- colMeans(cutof)

################################################################################
## Generate the simulated data
################################################################################
for (i in 1:n){
  arg1 = list(mu=mu1[i], sigma2=Sigma1, shape=shape1, nu=nu)
  arg2 = list(mu=mu2[i], sigma2=Sigma2, shape=shape2, nu=nu)
  compl= rmix(1, pii, "Skew.t", list(arg1,arg2), cluster=TRUE) # change to ("t", "Skew.t", "Skew.cn", "Skew.slash",
  #"Skew.normal", "Normal" from the mixsmsns package
  y[i] = compl$y
  gru[i] = compl$cluster
}

aa=sort(y,decreasing=FALSE)
cc=matrix(1,n,1)*(y<=cutoff)
y[cc==1]=cutoff

## creating censored observations    (right censoring)
# aa=sort(y,decreasing=TRUE)
# cc=matrix(1,n,1)*(y>=cutoff)
# y[cc==1]=cutoff

################################################################################
# Fitting the data using the ECME algorithm
# FM-NCR, FM-TCR, FM-SNCR, FM-STCR
################################################################################

#normal
fitN<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = 2, get.init = TRUE, criteria = TRUE, group = TRUE, family = "Normal", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)

# Student -t
fitT<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=4, g = 2, get.init = TRUE, criteria = TRUE, group = TRUE, family = "T", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)

# Skew normal
fitSN<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = 2, get.init = TRUE, criteria = TRUE, group = TRUE, family = "SN", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)

# Skew T      
fitST<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=4, g = 2, get.init = TRUE, criteria = TRUE, group = TRUE, family = "ST", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)



