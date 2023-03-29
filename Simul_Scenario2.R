
rm(list=ls(all=TRUE))

source("Moment_SMSNT.R")
source("Integral_nu_float.R")
source("algEM.fmr.ST.crFinal.R")

#install.packages("mixsmsn")
library(mixsmsn)

################################################################################
# Simulation Studies
# Simulation 2
# Three different configuration regression has considered - Parallel, Concurrent, and Trigonometric
################################################################################
# Case1: Parallel 
################################################################################

n<-1000 # n : 100, 200, 500, 1000

x<-cbind(1, runif(n,-1,3))

y<-matrix(0,n,1)
gru<-matrix(0,n,1)

betas1<-c(0, 1)
Sigma1 <- 0.1
mu1 <- x%*%betas1
shape1<- 5

betas2<-c(3,1)
Sigma2 <- 0.2
mu2 <- x%*%betas2
shape2<- -3

betas3<-c(-3,1)
Sigma3 <- 0.1
mu3 <- x%*%betas3
shape3<- -1

nu<-1.1
pii<-c(0.3, 0.5, 0.2)
################################################################################
# Case2: Concurrent
################################################################################

x<-cbind(1, runif(n,-1,3))

y<-matrix(0,n,1)
gru<-matrix(0,n,1)

betas1<-c(-1, 1)
Sigma1 <- 0.1 
mu1 <- x%*%betas1
shape1<- 5

betas2<-c(3,-1)
Sigma2 <- 0.2
mu2 <- x%*%betas2
shape2<- -3

betas3<-c(3,1)
Sigma3 <- 0.1
mu3 <- x%*%betas3
shape3<- -1

nu<-1.1
pii<-c(0.3, 0.5, 0.2)
perc<-0.05 # Censoring level = 0.05, 0.1, 0.2
################################################################################
# Case3: Trigonometric
################################################################################

x<-runif(n,0,1)
x<-cbind(1, x, sin(6*pi*x))

y<-matrix(0,n,1)
gru<-matrix(0,n,1)

betas1<-c(-3, 1, 1)
Sigma1 <- 0.1
mu1 <- x%*%betas1
shape1<- 5

betas2<-c(3,1, 1)
Sigma2 <- 0.2
mu2 <- x%*%betas2
shape2<- -3

betas3<-c(0, 1, 1)
Sigma3 <- 0.1
mu3 <- x%*%betas3
shape3<- -1

nu<-1.1
pii<-c(0.3, 0.5, 0.2)
################################################################################
## Generate the simulated data
################################################################################
for (i in 1:n){
  arg1 = list(mu=mu1[i], sigma2=Sigma1, shape=shape1, nu=nu)
  arg2 = list(mu=mu2[i], sigma2=Sigma2, shape=shape2, nu=nu)
  arg3 = list(mu=mu3[i], sigma2=Sigma3, shape=shape3, nu=nu)
  #compl= rmix(1, pii, "Skew.t", list(arg1,arg2), cluster=TRUE) 
  compl= rmix(1, pii, "Skew.slash", list(arg1,arg2,arg3), cluster=TRUE) # change to ("t", "Skew.t", "Skew.cn", "Skew.slash",
  #"Skew.normal", "Normal" from the mixsmsns package
  y[i] = compl$y
  gru[i] = compl$cluster
}

## creating censored observations (left censoring)
aa=sort(y,decreasing=FALSE)
cutof<-aa[ceiling(perc*n)]
cc=matrix(1,n,1)*(y<=cutof)
y[cc==1]=cutof

## creating censored observations (right censoring)
# aa=sort(y,decreasing=FALSE)
# cutof<-aa[ceiling(perc[k]*n)]
# cc=matrix(1,n,1)*(y>=cutof)
# y[cc==1]=cutof

################################################################################
# Fitting the data using the ECME algorithm
# FM-NCR, FM-TCR, FM-SNCR, FM-STCR
################################################################################

#normal
fitN<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "Normal", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)

# Student -t
fitT<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=4, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "T", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)

# Skew normal
fitSN<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "SN", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)

# Skew T      
fitST<-EM.skewCens.mixR(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=4, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "ST", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL,aitken = TRUE, IM=TRUE)



