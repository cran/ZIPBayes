ZIPMErr <- function(Ystar, Covarmainphi, Covarmainmu, Covarplus, Covarminus,
                                betaphi, betamu, alphaplus, alphaminus,
                                Uibound = c(7,11),
                                priorgamma, priormu, priorSigma, propsigmaphi,  propsigmamu = propsigmaphi, propsigmaplus = propsigmaphi,  propsigmaminus = propsigmaphi, 
                                seed = 1, nmcmc = 500){
  set.seed(seed)
  nsample <- length(Ystar)
  
  Zplus <- rep(0,nsample)
  Zminus<- rep(0,nsample)
  
  mu1 <-  as.matrix(Covarmainphi) %*% t(t(betaphi))
  # phi <- cexpexp(as.matrix(Covarmainphi) %*% t(t(betaphi)))
  mu2 <- as.matrix(Covarmainmu) %*% t(t(betamu))
  muZplus <- as.matrix(Covarplus) %*% t(t(alphaplus)) 
  muZminus <- as.matrix(Covarminus) %*% t(t(alphaminus))
  U1 <- rep(0,nsample)
  U2 <- rep(0,nsample)
  
  Y <- Ystar - Zplus + Zminus
  
  beta1record <- betaphi
  beta2record <- betamu
  alphaplusrecord <- alphaplus
  alphaminusrecord <- alphaminus
  for (repk in 1:nmcmc){
    
    ## Adjust the proposal standard error 
    if (repk %% 500 ==0){
      cat("MCMC interations:" , repk,"\n")
    }
    
    ### Step 1: Measurement Error: Positive Process
    for (t in 1:nsample){
      NewGen <-  ZI_GenerateBigJoint(Ystar[t], Uibound[1], Uibound[2],  mu1[t], mu2[t], muZplus[t],  muZminus[t])
      Zplus[t] <- NewGen[1]
      Zminus[t] <- NewGen[2]
      Y[t] <- NewGen[3]
      U1[t] <- NewGen[4]
      U2[t] <- NewGen[5]
    }
    
    alphaplussafe <- alphaplus
    
    alphaplus <- ZI_GenerateBetaMetro(alphaplus, Zplus, as.matrix(Covarplus), propsigmaplus,  priorgamma)
    
    muZplus <- as.matrix(Covarplus) %*% t(t(alphaplus))
    
    
    # ### Step 2: Measurement Error: Negative Process
    
    alphaminus <- ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminus),  propsigma= propsigmaminus,  priormu,  priorSigmas=priorSigma)
    
    muZminus <- as.matrix(Covarminus) %*% t(t(alphaminus))
    
    # ### Step 3: Main Model
    
    ### Step 3.1 Update betaphi and beta 2
    
    betaphi <- ZI_GenerateBetaMetro(Par=betaphi, U=U1, Covar=as.matrix(Covarmainphi), propsigma=propsigmaphi, priorgamma)
    mu1 <-  as.matrix(Covarmainphi) %*% t(t(betaphi))
   
    betamu <- ZI_GenerateBetaMetro(Par=betamu, U=U2, Covar=as.matrix(Covarmainmu), propsigma=propsigmamu, priorgamma)
    mu2 <- as.matrix(Covarmainmu) %*% t(t(betamu))
    
    
    beta1record <- rbind(beta1record,betaphi)
    beta2record <- rbind(beta2record,betamu)
    alphaplusrecord <- rbind(alphaplusrecord,alphaplus)
    alphaminusrecord <- rbind(alphaminusrecord,alphaminus)
  }
  
  BayesResults <- list(betaphi_trace = beta1record,
                       betamu_trace = beta2record, 
                       alphaplus_trace = alphaplusrecord, 
                       alphaminus_trace = alphaminusrecord)
  
  class(BayesResults) <- "ZIPBayes"
  
  return(BayesResults)
}
