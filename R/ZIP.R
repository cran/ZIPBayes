ZIP <- function(Y, Covarmainphi, Covarmainmu, 
                         betaphi, betamu, 
                         priorgamma, propsigmaphi,  propsigmamu = propsigmaphi,
                         seed = 1, nmcmc = 500){
  set.seed(seed)
  nsample <- length(Y)
  
  
  mu1 <-  as.matrix(Covarmainphi) %*% t(t(betaphi))
  phi <- cexpexp(as.matrix(Covarmainphi) %*% t(t(betaphi)))
  mu2 <- as.matrix(Covarmainmu) %*% t(t(betamu))
  U2 <- rep(1,length(Y)) 
  U1 <- rep(1,length(Y)) 
  
  beta1record <- betaphi
  beta2record <- betamu
 
  for (repk in 1:nmcmc){
    if (repk %% 500 ==0){
      cat("MCMC interations:" , repk,"\n")
    }
    
    # ### Step 3: Main Model
    
    ### Step 3.1 Data Augmentation
    U1 <- ZI_GenerateU1(Y,U2,exp(mu1))
    U2 <- ZI_GenerateU2(Y,U1,exp(mu2))
    
    ### Step 3.2 Update betaphi and beta 2
    betaphi <- ZI_GenerateBetaMetro(Par=betaphi, U=U1, Covar=as.matrix(Covarmainphi), propsigma=propsigmaphi, priorgamma)
    mu1 <-  as.matrix(Covarmainphi) %*% t(t(betaphi))
    phi <- cexpexp(mu1)
    
    
    betamu <- ZI_GenerateBetaMetro(Par=betamu, U=U2, Covar=as.matrix(Covarmainmu), propsigma=propsigmamu, priorgamma)
    mu2 <- as.matrix(Covarmainmu) %*% t(t(betamu))
    
    
    beta1record <- rbind(beta1record,betaphi)
    beta2record <- rbind(beta2record,betamu)
  }
  
  BayesResults <- list(betaphi_trace = beta1record, betamu_trace = beta2record)
  
  # results <- apply(BayesResults, MARGIN =2, FUN = function(x){
  #   xpure <- purifyseq(x, 500, 1)
  #   mean <- mean(xpure,na.rm=T)
  #   sd <- sd(xpure,na.rm=T)
  #   median <- median(xpure,na.rm=T)
  #   CrI <- quantile(xpure, c(0.025,0.975))
  #   HDP <- getHDP(xpure,0.05)
  #   return(c(mean,median,sd,CrI,HDP))
  # })
  
  class(BayesResults) <- "ZIPBayes"
  
  return(BayesResults)
}
