#library(mvtnorm)

#location-scale student-T distribution pdf
dent<-function(x,mu,sigma2,nu){
  z<-(x-mu)/sqrt(sigma2)
  aa<-1/sqrt(sigma2)*dt(x = z,df = nu)
  if(length(which(aa == 0)) > 0) aa[which(aa == 0)] <- .Machine$double.xmin
  return(aa)
}

#########################################################################

#location-scale student-T distribution pcf
pent<-function(x,mu,sigma2,nu){
aa<-pt((x-mu)/sqrt(sigma2),nu)
if(length(which(aa == 0)) > 0) aa[which(aa == 0)] <- .Machine$double.xmin
  return(aa)
}

##########################################################################

acumt2 = function(a = NULL,b,mu,Sigma,nu){
  Sigma21 = c(Sigma[-1,-1] - Sigma[-1,1]%*%solve(Sigma[1,1])%*%Sigma[1,-1])
  expab = function(x,y){
    delta1 = c((x - mu[1])^2/Sigma[1,1])
    mu21 = mu[-1] + Sigma[-1,1]%*%solve(Sigma[1,1])%*%(x - mu[1])
    return(pent(y,mu21,Sigma21*(nu+delta1)/(nu+1),nu+1))
  }
  func = function(x,val){
    val*dent(x,mu[1],Sigma[1,1],nu)
  }  
  if(all(is.infinite(a)) | is.null(a)){
  f = function(x){func(x,expab(x,b[2]))}
    return(integrate(Vectorize(f),lower = -Inf,upper = b[1])$value)
    #return(integrate(f = function(x) func(x,expab(x,b[2])),lower = mu[1]-37*sqrt(Sigma[1,1]),upper = b[1])$value)
  }else{
    return(integrate(f = function(x) func(x,expab(x,b[2])),lower = a[1],upper = b[1])$value - 
             integrate(f = function(x) func(x,expab(x,a[2])),lower = a[1],upper = b[1])$value)
  }
}

##########################################################################
#TEST
##########################################################################

# a = c(-1,1)
# b = c(2,3)
# mu = as.matrix(c(1,2))
# Sigma = matrix(c(2,-0.5,-0.5,1),2,2)
# 
# nu = 3
# acumt2(a,b,mu,Sigma,nu)
# pmvt(lower = c(a-mu),upper = c(b-mu),df = nu,sigma = Sigma)[1]

# 
# nu = 3.27
# acumt2(a,b,mu,Sigma,nu)
# pmvt(lower = c(a-mu),upper = c(b-mu),df = nu,sigma = Sigma)[1]
# 
# #CDF
# a = c(-Inf,-Inf)
# acumt2(a,b,mu,Sigma,nu)
# #or a = NULL
# acumt2(b = b,mu = mu,Sigma = Sigma,nu = nu)



#a = c(-0.8,-0.7,-0.6)
#b = c(0.5,0.6,0.7)
#mu = c(0.1,0.2,0.3)
#Sigma = matrix(data = c(1,0.2,0.3,0.2,1,0.4,0.3,0.4,1),
#pmvnormt(lower = a,upper = b,mean = mu,sigma = Sigma,nu = nu)
