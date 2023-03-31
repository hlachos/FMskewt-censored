
rm(list=ls(all=TRUE))

source("Moment_SMSNT.R")
source("Integral_nu_float.R")
source("algEM.fmr.ST.crFinal.R")

#################################################################################
### Aplication  WAGE RATES DATA SET used in Zeller Et al. 2019 (ADAC)
#################################################################################
# install.packages("SMNCensReg")
library(SMNCensReg) # Containing Wage Rates Data

data(wage.rates)
y <- wage.rates$hours/1000

# Check the density of the response variable
lattice:: densityplot(y)
hist(y, prob=TRUE, main="")
lines(density(y))

x <- cbind(1,wage.rates$educ,wage.rates$age,wage.rates$exper,wage.rates$expersq)

#censure
cc <- c(rep(0,428),rep(1,325))
n <- length(y)

# Normal
fitN<-EM.skewCens.mixR(cc, -y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = 2, get.init = TRUE, criteria = TRUE, group = FALSE, family = "Normal", error = 0.000001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL, aitken = TRUE, IM=TRUE)

# Student -t 
fitT<-EM.skewCens.mixR(cc, -y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=3, g = 2, get.init = TRUE, criteria = TRUE, group = FALSE, family = "T", error = 0.000001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL, aitken = TRUE, IM=TRUE)

# Skew normal 
fitSN<-EM.skewCens.mixR(cc, -y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = 2, get.init = TRUE, criteria = TRUE, group = FALSE, family = "SN", error =  0.000001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL, aitken = TRUE, IM=TRUE)

# Skew T 
fitST<-EM.skewCens.mixR(cc, -y, x, Abetas = fitSN$Abetas, sigma2 = fitSN$sigma2, shape = fitSN$shape, pii = fitSN$pii, nu=3,  g = 2, get.init = FALSE, criteria = TRUE, group = FALSE, family = "ST", error = 0.00001, iter.max = 1000, obs.prob= FALSE, kmeans.param = NULL, aitken = TRUE, IM=FALSE)
plot(fitST$loglikT, xlab="Iteration", ylab="Log-likelihood", main="ST", type="l")

