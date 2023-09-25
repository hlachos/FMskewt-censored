 ####################################################################################
 ### EM algorithm for Finite Mixtures of Censored Data under ST, SN, T and Normal
 ################################################################################
 
 EM.skewCens.mixR <- function(cc, y, x, Abetas = NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL, g = NULL, get.init = TRUE, criteria = TRUE, group = FALSE, family = "Normal", error = 0.00001, iter.max = 100, obs.prob= FALSE, kmeans.param = NULL, aitken = TRUE, IM=TRUE){
  #y:  data of size  n
  ## Abetas: regresion parameters
  #, sigma2,  pii: Initial values for the EM algorithm. 
  #nu: valor inicial para o nu (no caso da CN deve ser um vetor bidimensional com valores entre 0 e 1)
  # get.init = TRUE initial values are computed in the EM 
  #cluster: TRUE ou FALSE - Which component belong each observation
  #type: c("ST","T","SN","Normal") - T: fit a mixture of Student-T
  #                                       - ST: fit a mixture of ST
  #                                       - SN: fit a mixture ofSN
  #                                       - Normal: it a mixture of Normals
  #error: stopping criteria
  if(ncol(as.matrix(y)) > 1) stop("This function is only for univariate response y!") 
  if((family != "ST") && (family != "Normal") && (family != "T") && (family != "SN")) stop(paste("Family",family,"not recognized.\n", sep=" "))
  if((length(g) == 0) && ((length(sigma2)==0) || (length(shape)==0) || (length(pii)==0) || (ncol(Abetas)==0) ))  stop("The model is not specified correctly.\n")
  if(get.init == FALSE){
       if((length(pii) != length(sigma2)) || (ncol(Abetas) != length(pii))) stop("The size of the initial values are not compatibles.\n")
       if((family == "ST" || family == "Normal" || family == "T" || family == "SN") & (length(pii) != ncol(Abetas))) stop("The size of the initial values are not compatibles.\n")
       if(sum(pii) != 1) stop("probability of pii does not sum to 1.\n")
  }  
  
  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  p<-ncol(x)
  n <- length(y)


  if (get.init == TRUE){
  if(length(g) == 0) stop("g is not specified correctly.\n")

    k.iter.max <- 50
    n.start <- 1
    algorithm <- "Hartigan-Wong"
    if(length(kmeans.param) > 0){
       if(length(kmeans.param$iter.max) > 0 ) k.iter.max <- kmeans.param$iter.max
       if(length(kmeans.param$n.start) > 0 ) n.start <- kmeans.param$n.start
       if(length(kmeans.param$algorithm) > 0 ) algorithm <- kmeans.param$algorithm
    }

    if(g > 1){
      nu<-3
      init <- kmeans(y,g,k.iter.max,n.start,algorithm)
      pii <- init$size/length(y)
      sigma2 <- c()
      Abetas<-matrix(0,p,g)
      shape <- c()
     for (j in 1:g){
     dd <- init$cluster 
     Abetas[,j] <- solve(t(x[dd==j,])%*%x[dd==j,])%*%t(x[dd==j,])%*%y[dd==j]
      yr<- y[dd==j]-x[dd==j,]%*%Abetas[,j]
      sigma2[j]<-sum(yr^2)/init$size[j]
      m3 <- (1/init$size[j])*sum( (y[init$cluster == j] - x[dd==j,]%*%Abetas[,j])^3 )
      shape[j] <- 3*sign(m3/sigma2[j]^(3/2))
      }
      
    }
    
   else{
     Abetas<-solve(t(x)%*%x)%*%t(x)%*%y
     yr<-y-x%*%Abetas
     sigma2 <- var(yr)
     m3 <- (1/length(yr))*sum( (yr)^3 )
     shape <- 2*sign(m3/sigma2^(3/2))
     pii <- 1
     nu<-4
   }
  }

  #Lim1 <- y
  #Lim2 <- rep(Inf,n)
   loglikT<-numeric()
     
################################################################################
###                                     Normal
################################################################################
if (family == "Normal"){

      mu<- matrix(0,n,g)
      delta <- matrix(0,g,1)
      Delta <- matrix(0,g,1)
      Gama <- matrix(0,g,1)
      for (k in 1:g){

        mu[,k]<- x%*%Abetas[,k]
       shape[k] <-0
      delta[k] <- 0#shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- 0# sqrt(sigma2[k])*delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2   
     }
      
      criterio <- 1
      count <- 0

      lk = lk1 = lk2  <- sum(log(d.mixedN(cc, y, pii, mu, sigma2)))
      
      
      
      while((criterio > error) && (count <= iter.max)){
      
        count <- count + 1
         print(count)
   
        tal <- matrix(0, n, g)
        #soma1<-matrix(0, p,1)
        #soma2<-matrix(0, p, p)
        

        for (k in 1:g){
          ### E-step: calculando os momentos
          
     dj <- ((y - mu[, k])/sqrt(sigma2[k]))^2
      Mtij2 <- 1/(1 + (Delta[k]^2)*(Gama[k]^(-1)))
      Mtij <- sqrt(Mtij2)
      mutij <- Mtij2*Delta[k]*(Gama[k]^(-1))*(y - mu[, k])
      A <- mutij / Mtij
      cnu<-2/sqrt(2*pi*(1+shape[k]^2))
      prob <- pnorm(mutij/Mtij)
      if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
      E= dnorm(mutij/Mtij)/prob
      u= rep(1, n)
      
      S1 <- u
      S2 <- (mutij*u + Mtij*E)
      S3 <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00<-S1
      E01<-y*S1
      E02<-y^2*S1
      E10<-S2
      E20<-S3
      E11<-y*S2
      sigma2s<- (1/(1+shape[k]^2))*sigma2[k]
      
      Aux11<-cdfSNI(y, mu[, k], sigma2s, 0, NULL, type = "Normal")
      Aux2<- cdfSNI(y, mu[, k], sigma2[k], 0, NULL, type = "Normal")
        
      mu1<-mu[cc==1,k]
      y1<-y[cc==1]
      np<-length(mu1)
      aux1MomW<-aux2MomS<-matrix(0,np,2)
       
      for(j in 1:np){
        A1a<-meanvarTMD( -Inf,y1[j],mu1[j],sigma2[k],dist = "normal")
        #MomenSNI(mu1[j], sigma2[k], 0, NULL, delta=NULL,Lim1=-Inf,Lim2=y1[j], type="SN") #CalMom(mu,sigma2,nu,y,type)
        A2a<- meanvarTMD( -Inf,y1[j],mu1[j],sigma2s,dist = "normal")
        #MomenSNI(mu1[j], sigma2s, 0, NULL, delta=NULL, Lim1=-Inf, Lim2=y1[j], type="SN")
        aux1MomW[j,]<-c(A1a$mean,A1a$EYY)
        aux2MomS[j,]<-c(A2a$mean,A2a$EYY)
                }
      
      Wphi<- Aux11/Aux2
      
      E00Aux<-E01Aux<-E02Aux<-WW1<-u
      E01Aux[cc==1]<-aux1MomW[,1]
      E02Aux[cc==1]<-aux1MomW[,2]
      WW1[cc==1]<-aux2MomS[,1]
      
      
      E10Aux<- Delta[k]/(Gama[k] + Delta[k]^2)*(E01Aux-E00Aux*mu[, k])+Mtij*cnu*Wphi
      E20Aux<- (Delta[k]/(Gama[k] + Delta[k]^2))^2*(E02Aux-2*E01Aux*mu[, k]+mu[, k]^2*E00Aux)+Mtij*(Delta[k]/(Gama[k] + Delta[k]^2))*cnu*Wphi*(WW1-mu[, k])+ Mtij2
      E11Aux<-Delta[k]/(Gama[k] + Delta[k]^2)*(E02Aux-E01Aux*mu[, k])+ Mtij*Wphi*cnu*WW1
      
      
      E00[cc==1]<- E00Aux[cc==1]
      E01[cc==1]<- E01Aux[cc==1]
      E02[cc==1]<- E02Aux[cc==1]
      E10[cc==1]<- E10Aux[cc==1]
      E20[cc==1]<-E20Aux[cc==1]
      E11[cc==1]<-E11Aux[cc==1]
      
        
                          

          d1 <- dNormal(cc, y, mu[,k], sigma2[k])
          
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedN(cc, y, pii, mu, sigma2)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
      
          ### M-step: atualizar os parametros ###
          
          pii[k] <- (1/n)*sum(tal[,k])
          
                   
      Abetas[,k]<- solve(t(x)%*%diag(tal[,k]*E00)%*%x)%*%(t(x)%*%(tal[,k]*E01))
      mu[,k]<-x%*%Abetas[,k]
      Delta[k]<-0#sum(E11-E10*mu)/sum(E20)
      Gama[k]<-sum(tal[,k]*E02-2*tal[,k]*E01*mu[, k]+tal[,k]*E00*mu[, k]^2)/sum(tal[,k])

      sigma2[k] <- Gama[k] + Delta[k]^2
      shape[k] <- 0#((sigma2[k]^(-1/2))*Delta[k])/(sqrt(1 - (Delta[k]^2)*(sigma2[k]^(-1))))
      
       #print(c(Abetas,sigma2,shape))    
        }
        

        
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }
        
        
      if (pii[1]< 0.5 & g==2){
      
      mu<-cbind(mu[,2],mu[,1])
       Abetas<-cbind(Abetas[,2], Abetas[,1])
       pii<- as.vector(c(pii[2], pii[1]))
       sigma2<-as.vector(c(sigma2[2], sigma2[1]))
      shape<-as.vector(c(shape[2], shape[1]))
      }
      

         

        lk <- sum(log(d.mixedN(cc, y, pii, mu, sigma2) ))
        #criterio <- abs(lk1/lk-1)
        #lk<-lk1
        loglikT[count]<-lk
                 if(aitken==TRUE)
      {
        lk3         <- lk

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)#;print(tmp)
          criterio  <- abs(tmp2 - lk3)
        }
        lk2         <- lk3
      }else{
        lk1        <- lk
        criterio   <- abs(lk1/lk2-1)
        lk2        <- lk1
      }
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
   #    icl <- 0
   #    for (j in 1:g) icl<-icl+sum(log(pii[j]*dNormal(cc, y, mu[,j], sigma2[j])))
    }
    if(IM==TRUE){
    for(k in 1:g){
    d1 <- dNormal(cc, y, mu[,k], sigma2[k])
          
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedN(cc, y, pii, mu, sigma2)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
      
    }
  desvios <-im.fmr.st.cr (cc, y, x, Abetas, sigma2, shape, pii, NULL,tal=tal, "Normal")$EP
             }
             else desvios=0
}

################################################################################
###                                     Student-t
################################################################################

if (family == "T"){

      
      mu<- matrix(0,n,g)
      delta <- matrix(0,g,1)
      Delta <- matrix(0,g,1)
      Gama <- matrix(0,g,1)
      
      for (k in 1:g){

       mu[,k]<- x%*%Abetas[,k]
      shape[k]<-0  
      delta[k] <- 0
      Delta[k] <- 0
      Gama[k] <- sigma2[k] - Delta[k]^2
      
      
     }
      
      criterio <- 1
      count <- 0

      lk = lk1 = lk2  <- sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu)))
      
      while((criterio > error) && (count <= iter.max)){
      
        count <- count + 1
         print(count)
   
        tal <- matrix(0, n, g)
        
        
           
        for (k in 1:g){
          ### E-step: calculando os momentos
          
        dj <- ((y - mu[, k])/sqrt(sigma2[k]))^2
      Mtij2 <- 1/(1 + (Delta[k]^2)*(Gama[k]^(-1)))
      Mtij <- sqrt(Mtij2)
      mutij <- Mtij2*Delta[k]*(Gama[k]^(-1))*(y - mu[, k])
      A <- mutij / Mtij
      cnu<-2*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi*(1+shape[k]^2)))
      E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2[k])*dt.ls(y, mu[, k], sigma2[k],shape[k] ,nu))
      u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2[k])*dt.ls(y, mu[, k], sigma2[k],shape[k] ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)
      
      S1 <- u
      S2 <- (mutij*u + Mtij*E)
      S3 <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00<-S1
      E01<-y*S1
      E02<-y^2*S1
      E10<-S2
      E20<-S3
      E11<-y*S2
      
      sigma2s<- nu/(nu+2)*sigma2[k]
      sigma2ss<- nu/((nu+1)*(1+shape[k]^2))*sigma2[k]
      
      Aux1<- cdfSNI(y, mu[, k], sigma2s, 0, nu+2, type = "T") 
      Aux11<-cdfSNI(y, mu[, k], sigma2ss, 0, nu+1, type = "T")
      Aux2<- cdfSNI(y, mu[, k], sigma2[k], 0, nu, type = "T")
      
      mu1<-mu[cc==1,k]
      y1<-y[cc==1]
      np<-length(mu1)
      aux1MomW<-aux2MomS<-matrix(0,np,2)
      
      # aux1MomW<-aux2MomS<-matrix(0,n,2)
      
      for(j in 1:np){      
        A1a<-   meanvarTMD(-Inf,y1[j],mu1[j],sigma2s,nu=nu+2,dist = "t")
        #MomenSNI(mu1[j],sigma2s,0,nu+2,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="ST") #CalMom(mu,sigma2,nu,y,type)
        A2a<-  meanvarTMD( -Inf,y1[j],mu1[j],sigma2ss,nu=nu+1,dist = "t")
        #MomenSNI(mu1[j],sigma2ss,0,nu+1,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="ST")
        aux1MomW[j,]<-c(A1a$mean,A1a$EYY)
        aux2MomS[j,]<-c(A2a$mean,A2a$EYY)
        
      }
      
      P1aux<-P2aux<-WW<-u
      P1aux[cc==1]<-aux1MomW[,1]
      P2aux[cc==1]<-aux1MomW[,2]
      WW[cc==1]<-aux2MomS[,1]
      
      Wphi<- Aux11/Aux2
      E00Aux<-Aux1/Aux2
      E01Aux<-Aux1/Aux2*P1aux
      E02Aux<-Aux1/Aux2* P2aux
      
      E10Aux<- Delta[k]/(Gama[k] + Delta[k]^2)*(E01Aux-E00Aux*mu[, k])+Mtij*cnu*Wphi
      E20Aux<- (Delta[k]/(Gama[k] + Delta[k]^2))^2*(E02Aux-2*E01Aux*mu[, k]+mu[, k]^2*E00Aux)+Mtij*(Delta[k]/(Gama[k] + Delta[k]^2))*cnu*Wphi*(WW-mu[, k])+ Mtij2
      E11Aux<-Delta[k]/(Gama[k] + Delta[k]^2)*(E02Aux-E01Aux*mu[, k])+ Mtij*Wphi*cnu*WW
      
      
      
      E00[cc==1]<- E00Aux[cc==1]
      E01[cc==1]<- E01Aux[cc==1]
      E02[cc==1]<- E02Aux[cc==1]
      E10[cc==1]<- E10Aux[cc==1]
      E20[cc==1]<-E20Aux[cc==1]
      E11[cc==1]<-E11Aux[cc==1]
      
      
          
          
          d1 <- dT(cc, y, mu[,k], sigma2[k],nu)
          
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedT(cc, y, pii, mu, sigma2,nu)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
      
          ### M-step: atualizar os parametros ###
          
          pii[k] <- (1/n)*sum(tal[,k])
          
      Abetas[,k]<- solve(t(x)%*%diag(c(E00*tal[,k]))%*%x)%*%(t(x)%*%matrix(E01*tal[,k],n,1))
      mu[,k]<-x%*%Abetas[,k]
      Delta[k]<-0#sum(E11-E10*mu)/sum(E20)
      Gama[k]<-sum(E02*tal[,k]-2*E01*tal[,k]*mu[, k]+E00*tal[,k]*mu[, k]^2)/sum(tal[,k])		
      sigma2[k] <- Gama[k] + Delta[k]^2
      shape[k] <- 0#((sigma2[k]^(-1/2))*Delta[k])/(sqrt(1 - (Delta[k]^2)*(sigma2[k]^(-1))))
   
   
   
        }
        

        
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }
        
        
        if (pii[1]< 0.5 & g==2){
      
       mu<-cbind(mu[,2],mu[,1])
       Abetas<-cbind(Abetas[,2], Abetas[,1])
       pii<- as.vector(c(pii[2], pii[1]))
       sigma2<-as.vector(c(sigma2[2], sigma2[1]))
       shape<-as.vector(c(shape[2], shape[1]))
       
      }
    

         
        ft <- function(nu)sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu)))
		#	  nu <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum 
       nu<-optimize(ft, c(2.01,150), tol = 0.000001, maximum = TRUE)$maximum  
        lk <- ft(nu)
        loglikT[count]<-lk
       # criterio <- abs(lk1/lk-1)
       # lk<-lk1

               if(aitken==TRUE)
      {
        lk3         <- lk

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)#;print(tmp)
          criterio  <- abs(tmp2 - lk3)
        }
        lk2         <- lk3
      }else{
        lk1        <- lk
        criterio   <- abs(lk1/lk2-1)
        lk2        <- lk1
      }
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
    #   icl <- 0
    #   for (j in 1:g) icl<-icl+sum(log(pii[j]*dT(cc, y, mu[,j], sigma2[j],nu)))
    }
    if(IM==TRUE){
              for(k in 1:g){
          d1 <- dT(cc, y, mu[,k], sigma2[k],nu)
          
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedT(cc, y, pii, mu, sigma2,nu)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
          
          }
      
  desvios  <-im.fmr.st.cr (cc, y, x, Abetas, sigma2, shape, pii, nu, tal=tal, "T")$EP
             }
             else desvios=0
 
}


################################################################################
###                                     Skew-T
################################################################################

if (family == "ST"){

     
       mu<- matrix(0,n,g)
      delta <- matrix(0,g,1)
      Delta <- matrix(0,g,1)
      Gama <- matrix(0,g,1)
      
      for (k in 1:g){

      mu[,k]<- x%*%Abetas[,k]
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k])*delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      
     }
      
      criterio <- 1
      count <- 0

      lk = lk1 = lk2 <- sum(log(d.mixedSTC1(cc, y, pii, mu, sigma2, shape, nu)))
      
      while((criterio > error) && (count <= iter.max)){
      
        count <- count + 1
         print(count)
   
        tal <- matrix(0, n, g)
           
        for (k in 1:g){
          ### E-step: calculando os momentos
      
      
       dj <- ((y - mu[, k])/sqrt(sigma2[k]))^2
      Mtij2 <- 1/(1 + (Delta[k]^2)*(Gama[k]^(-1)))
      Mtij <- sqrt(Mtij2)
      mutij <- Mtij2*Delta[k]*(Gama[k]^(-1))*(y - mu[, k])
      A <- mutij / Mtij
      cnu<-2*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi*(1+shape[k]^2)))
      E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2[k])*dt.ls(y, mu[, k], sigma2[k],shape[k] ,nu))
      u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2[k])*dt.ls(y, mu[, k], sigma2[k],shape[k] ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)
      
      S1 <- u
      S2 <- (mutij*u + Mtij*E)
      S3 <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00<-S1
      E01<-y*S1
      E02<-y^2*S1
      E10<-S2
      E20<-S3
      E11<-y*S2
      ychap<-y
      sigma2s<- nu/(nu+2)*sigma2[k]
      sigma2ss<- nu/((nu+1)*(1+shape[k]^2))*sigma2[k]
      
      Aux1<- cdfSNI(y, mu[, k], sigma2s, shape[k], nu+2, type = "ST") 
      Aux11<-cdfSNI(y, mu[, k], sigma2ss, 0, nu+1, type = "ST")
      Aux2<- cdfSNI(y, mu[, k], sigma2[k], shape[k], nu, type = "ST")
      
      
      mu1<-mu[cc==1, k]
      y1<-y[cc==1]
      np<-length(mu1)
      aux1MomW<-aux2MomS<-aux3MomS<-matrix(0,np,2)
      
      
      for(j in 1:np){      
        A1a<-   meanvarTMD( -Inf,y1[j],mu1[j],sigma2s,shape[k], nu=nu+2,dist = "ST")
        #MomenSNI(mu1[j],sigma2s,shape[k],nu+2,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="ST") #CalMom(mu,sigma2,nu,y,type)
        A2a<-  meanvarTMD( -Inf,y1[j],mu1[j],sigma2ss,nu=nu+1,dist = "t")
        #MomenSNI(mu1[j],sigma2ss,0,nu+1,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="ST")
        A3a<-    meanvarTMD(-Inf,y1[j],mu1[j],sigma2[k],nu=nu,dist = "t")
        #MomenSNI(mu1[j],sigma2[k],0,nu,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="ST")
        aux1MomW[j,]<-c(A1a$mean,A1a$EYY)
        aux2MomS[j,]<-c(A2a$mean,A2a$EYY)
        aux3MomS[j,]<-c(A3a$mean,A3a$EYY)
        
        
      }
      
      P1aux<-P2aux<-WW<-u
      P1aux[cc==1]<-aux1MomW[,1]
      P2aux[cc==1]<-aux1MomW[,2]
      WW[cc==1]<-aux2MomS[,1]
      
      
      Wphi<- Aux11/Aux2
      E00Aux<-Aux1/Aux2
      E01Aux<-Aux1/Aux2*P1aux
      E02Aux<-Aux1/Aux2*P2aux
      
      E10Aux<- Delta[k]/(Gama[k] + Delta[k]^2)*(E01Aux-E00Aux*mu[, k])+Mtij*cnu*Wphi
      E20Aux<- (Delta[k]/(Gama[k] + Delta[k]^2))^2*(E02Aux-2*E01Aux*mu[, k]+mu[, k]^2*E00Aux)+Mtij*(Delta[k]/(Gama[k] + Delta[k]^2))*cnu*Wphi*(WW-mu[, k])+ Mtij2
      E11Aux<- Delta[k]/(Gama[k] + Delta[k]^2)*(E02Aux-E01Aux*mu[, k])+ Mtij*Wphi*cnu*WW
      
      ychap[cc==1]<-aux3MomS[,1]
      E00[cc==1]<- E00Aux[cc==1]        ## E[U]
      E01[cc==1]<- E01Aux[cc==1]        ## E[UY]  
      E02[cc==1]<- E02Aux[cc==1]        ## E[UY2]   
      E10[cc==1]<- E10Aux[cc==1]        ## E[UT]
      E20[cc==1]<-E20Aux[cc==1]         ## E[UT2]
      E11[cc==1]<-E11Aux[cc==1]         ## E[UTY]  
                 
          d1 <- dSTC1(cc, y, mu[,k], sigma2[k] ,shape[k], nu)
                    
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedSTC1(cc, y, pii, mu, sigma2, shape, nu)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
      
          ### M-step: atualizar os parametros ###
          
          pii[k] <- (1/n)*sum(tal[,k])
           Abetas[,k]<- solve(t(x)%*%diag(c(E00*tal[,k]))%*%x)%*%(t(x)%*%matrix(E01*tal[,k]-E10*tal[,k]*Delta[k],n,1))
      mu[,k]<-x%*%Abetas[,k]
      Delta[k]<-sum(E11*tal[,k]-tal[,k]*E10*mu[, k])/sum(E20*tal[,k])
      Gama[k]<-sum(E02*tal[,k]-2*E01*tal[,k]*mu[, k]+E00*tal[,k]*mu[, k]^2+Delta[k]^2*E20*tal[,k]-2*Delta[k]*E11*tal[,k]+2*Delta[k]*E10*tal[,k]*mu[, k])/sum(tal[,k])		
      sigma2[k] <- Gama[k] + Delta[k]^2
      shape[k] <- ((sigma2[k]^(-1/2))*Delta[k])/(sqrt(1 - (Delta[k]^2)*(sigma2[k]^(-1))))
       # print(c(Abetas,sigma2,shape,nu))   
        }
        

        
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }
        
        
        if (pii[1]< 0.5 & g==2){
      
        mu<-cbind(mu[,2],mu[,1])
       Abetas<-cbind(Abetas[,2], Abetas[,1])
       pii<- as.vector(c(pii[2], pii[1]))
       sigma2<-as.vector(c(sigma2[2], sigma2[1]))
       shape<-as.vector(c(shape[2], shape[1]))
       
      }
    

         
   ft <- function(nu)sum(log(d.mixedSTC1(cc, y, pii, mu, sigma2, shape, nu)))
			  #nu <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),lower = MIN_NU, upper=MAX_NU,maximum=TRUE,tol=TOLERANCIA)$maximum 
			  #nu<-optimize(ft, c(2.1,40), tol = 0.00001, maximum = TRUE)$maximum
                             
        #print(nu)
        #lk1 <- ft(nu)
       
       lk<-ft(nu)
         loglikT[count]<-lk
       # criterio <- abs(lk1/lk-1)
       # lk<-lk1
             if(aitken==TRUE)
      {
        lk3         <-lk

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)#;print(tmp)
          criterio  <- abs(tmp2 - lk3)
        }
        lk2         <- lk3
      }else{
        lk1        <- lk
        criterio   <- abs(lk1/lk2-1)
        lk2        <- lk1
      }
      }

    if (criteria == TRUE){
      cl <- apply(tal, 1, which.max)
   }

if(IM==TRUE){
for(k in 1:g){
         d1 <- dSTC1(cc, y, mu[,k], sigma2[k] ,shape[k], nu)
                    
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedSTC1(cc, y, pii, mu, sigma2, shape, nu)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
}

  desvios <-im.fmr.st.cr (cc, y, x, Abetas, sigma2, shape, pii, nu, tal=tal, "ST")$EP
             }
             else desvios=0
 
}

################################################################################
###                                     Skew Normal
################################################################################

if (family == "SN"){
 
      mu<- matrix(0,n,g)
      delta <- matrix(0,g,1)
      Delta <- matrix(0,g,1)
      Gama <- matrix(0,g,1)
      
      for (k in 1:g){

       
        mu[,k]<- x%*%Abetas[,k]
        
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k])*delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      
     }
      
      criterio <- 1
      count <- 0

      lk = lk1 = lk2 <- sum(log(d.mixedSNC1(cc, y, pii, mu, sigma2, shape)))
      
      while((criterio > error) && (count <= iter.max)){
      
        count <- count + 1
         print(count)
   
        tal <- matrix(0, n, g)
        
        
        for (k in 1:g){
          ### E-step: calculando os momentos
         
      dj <- ((y - mu[, k])/sqrt(sigma2[k]))^2
      Mtij2 <- 1/(1 + (Delta[k]^2)*(Gama[k]^(-1)))
      Mtij <- sqrt(Mtij2)
      
      mutij <- Mtij2*Delta[k]*(Gama[k]^(-1))*(y - mu[, k])
      A <- mutij / Mtij   
      cnu<- 2/sqrt(2*pi*(1+shape[k]^2))
      prob <- pnorm(mutij/Mtij)
      if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
      E= dnorm(mutij/Mtij)/prob
      u= rep(1, n)
      
      S1 <- u
      S2 <- (mutij*u + Mtij*E)
      S3 <- (mutij^2*u + Mtij2 + Mtij*mutij*E)
      
      E00<-S1
      E01<-y*S1
      E02<-y^2*S1
      E10<-S2
      E20<-S3
      E11<-y*S2
      ychap<-y
      #sigma2s<- nu/(nu+2)*sigma2
      sigma2s<- (1/(1+shape[k]^2))*sigma2[k]
      
      Aux11<-cdfSNI(y, mu[, k], sigma2s, 0, NULL, type = "Normal")
      Aux2<- cdfSNI(y, mu[, k], sigma2[k], shape[k], NULL, type = "SN")
      
      mu0<- mu[, k]
      mu1<- mu0[cc==1]
      y1<-y[cc==1]
      np<-length(mu1)
      aux1MomW<-aux2MomS<-aux3MomS<-matrix(0,np,2)
       
      for(j in 1:np){
        A1a<-  meanvarTMD( -Inf,y1[j],mu1[j],sigma2[k],shape[k],dist = "SN")
        #MomenSNI(mu1[j],sigma2[k],shape[k],NULL,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="SN") #CalMom(mu,sigma2,nu,y,type)
        A2a<- meanvarTMD( -Inf,y1[j],mu1[j],sigma2s,dist = "normal")
        #MomenSNI(mu1[j],sigma2s,0,NULL,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="SN")
        A3a<- meanvarTMD(-Inf,y1[j],mu1[j],sigma2[k],dist = "normal")
        #MomenSNI(mu1[j],sigma2[k],0,NULL,delta=NULL,Lim1=-Inf,Lim2=y1[j],type="SN")
        aux1MomW[j,]<-c(A1a$mean,A1a$EYY)#c(A1a$EUY1,A1a$EUY2)
        aux2MomS[j,]<-c(A2a$mean,A2a$EYY)
        aux3MomS[j,]<-c(A3a$mean,A3a$EYY)
        }
      
      Wphi<- Aux11/Aux2
      
      E00Aux<-E01Aux<-E02Aux<-WW1<-u
      E01Aux[cc==1]<-aux1MomW[,1]
      E02Aux[cc==1]<-aux1MomW[,2]
      WW1[cc==1]<-aux2MomS[,1]
      
      
      E10Aux<- Delta[k]/(Gama[k] + Delta[k]^2)*(E01Aux-E00Aux*mu[, k])+Mtij*cnu*Wphi
      E20Aux<- (Delta[k]/(Gama[k] + Delta[k]^2))^2*(E02Aux-2*E01Aux*mu[, k]+mu[, k]^2*E00Aux)+Mtij*(Delta[k]/(Gama[k] + Delta[k]^2))*cnu*Wphi*(WW1-mu[, k])+ Mtij2
      E11Aux<-Delta[k]/(Gama[k] + Delta[k]^2)*(E02Aux-E01Aux*mu[, k])+ Mtij*Wphi*cnu*WW1
      
      ychap[cc==1]<-aux3MomS[,1]
      E00[cc==1]<- E00Aux[cc==1]
      E01[cc==1]<- E01Aux[cc==1]
      E02[cc==1]<- E02Aux[cc==1]
      E10[cc==1]<- E10Aux[cc==1]
      E20[cc==1]<- E20Aux[cc==1]
      E11[cc==1]<- E11Aux[cc==1]
      
          
          d1 <- dSNC1(cc, y, mu[, k], sigma2[k], shape[k])
          
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedSNC1(cc, y, pii, mu, sigma2, shape)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2
      
          ### M-step: update of the parameters ###
          
          pii[k] <- (1/n)*sum(tal[,k])
          
          
       Abetas[,k]<- solve(t(x)%*%diag(c(E00*tal[,k]))%*%x)%*%(t(x)%*%matrix(E01*tal[,k]-E10*tal[,k]*Delta[k],n,1))
      mu[,k]<-x%*%Abetas[,k]
      Gama[k]<-sum(E02*tal[,k]-2*E01*tal[,k]*mu[, k]+E00*tal[,k]*mu[, k]^2+Delta[k]^2*E20*tal[,k]-2*Delta[k]*E11*tal[,k]+2*Delta[k]*E10*tal[,k]*mu[, k])/sum(tal[,k])		
      Delta[k]<-sum(E11*tal[,k]-tal[,k]*E10*mu[, k])/sum(E20*tal[,k])
      sigma2[k] <- Gama[k] + Delta[k]^2
      shape[k] <- ((sigma2[k]^(-1/2))*Delta[k])/(sqrt(1 - (Delta[k]^2)*(sigma2[k]^(-1))))
      #   print(c(Abetas,sigma2,shape,nu))
        }
        

        
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }
        
        
        if (pii[1]< 0.5 & g==2){
       mu<-cbind(mu[,2],mu[,1])
       Abetas<-cbind(Abetas[,2], Abetas[,1])
       pii<- as.vector(c(pii[2], pii[1]))
       sigma2<-as.vector(c(sigma2[2], sigma2[1]))
       shape<-as.vector(c(shape[2], shape[1]))
      }
    
        	  
       # lk1 <- sum(log(d.mixedSNC1(cc, y, pii, mu, sigma2,shape) ))
        #criterio <- abs(lk1/lk-1)
        #lk<-lk1

        lk<-sum(log(d.mixedSNC1(cc, y, pii, mu, sigma2,shape) ))
        loglikT[count]<-lk            
        if(aitken==TRUE)
      {
        lk3         <- lk

        if(count<2){ criterio <- abs(lk2 - lk3)/abs(lk3)
        }else {
          c         <- (lk3 - lk2)/(lk2 - lk1)
          tmp2      <- lk2 + (lk3 - lk2)/(1-c)#;print(tmp)
          criterio  <- abs(tmp2 - lk3)
        }
        lk2         <- lk3
      }else{
        lk1        <- lk
        criterio   <- abs(lk1/lk2-1)
        lk2        <- lk1
      }
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
  }
      if(IM==TRUE){
      for(k in 1:g){
           d1 <- dSNC1(cc, y, mu[, k], sigma2[k], shape[k])
          
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          
          d2 <- d.mixedSNC1(cc, y, pii, mu, sigma2, shape)
          
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin
          
          tal[,k] <- d1*pii[k]/d2}
      
  desvios  <- im.fmr.st.cr (cc, y, x, Abetas, sigma2, shape, pii, NULL, tal=tal, "SN")$EP
             }
             else desvios=0
}



  if (criteria == TRUE){       
    if(family == "Normal") d <- g*p+ g + (g-1) #Abeta + sigma + pi
    if(family == "SN") d <- g*p+ 2*g + (g-1) #Abeta + sigma + lambda + pi
    if(family == "T") d <- g*p+ g + (g-1)+ 1 #Abeta + sigma +pi +nu
    if(family == "ST") d <- g*p+ 2*g + (g-1)+ 1 #Abeta + sigma +lambda+ pi +nu
    aic <- -2*lk + 2*d
    bic <- -2*lk + log(n)*d
    edc <- -2*lk + 0.2*sqrt(n)*d
    obj.out <- list(Abetas = Abetas, sigma2 = sigma2, shape=shape, pii = pii, sd=desvios,  nu=nu,  loglik=lk, loglikT=loglikT ,aic = aic, bic = bic, edc = edc,  iter = count, n = length(y), group = cl)#, im.sdev=NULL)
  }
 else obj.out <- list(Abetas = Abetas, sigma2 = sigma2, shape=shape, pii = pii, sd=desvios, nu=nu, loglik=lk,  iter = count, n = length(y), group = apply(tal, 1, which.max))#, im.sdev=NULL)
# if (group == FALSE) obj.out <- obj.out[-(length(obj.out)-1)]
 if (group == FALSE) obj.out <- obj.out[-(length(obj.out))]
 if (obs.prob == TRUE){
     nam <- c()
     for (i in 1:ncol(tal)) nam <- c(nam,paste("Group ",i,sep=""))
     if(ncol(tal) == 1) dimnames(tal)[[2]] <- list(nam)
     if(ncol(tal) > 1) dimnames(tal)[[2]] <- nam
     obj.out$obs.prob <- tal
     if((ncol(tal) - 1) > 1) obj.out$obs.prob[,ncol(tal)] <- 1 - rowSums(obj.out$obs.prob[,1:(ncol(tal)-1)])
     else obj.out$obs.prob[,ncol(tal)] <- 1 - obj.out$obs.prob[,1]
     obj.out$obs.prob[which(obj.out$obs.prob[,ncol(tal)] <0),ncol(tal)] <- 0.0000000000
     obj.out$obs.prob <- round(obj.out$obs.prob,10)
 }
 class(obj.out) <- family

 
# else obj.out <- obj.out[-length(obj.out)] 
 
 class(obj.out) <- family
 obj.out
}

