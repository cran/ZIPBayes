ZIPExt <- function(Ystar, Covarmainphi, Covarmainmu, Covarplus, Covarminus,
           Ystarval, Yval, Covarvalplus, Covarvalminus,
           betaphi, betamu, alphaplus, alphaminus,
           Uibound = c(7,11),
           priorgamma, priormu, priorSigma, 
           propsigmaphi,  propsigmamu = propsigmaphi, propsigmaplus = propsigmaphi,  propsigmaminus = propsigmaphi,
           seed = 1, nmcmc = 500){
  
    set.seed(seed)
    nmain <- length(Ystar)
    nval <- length(Ystarval)
    
    nsample <- nmain + nval
    
    Zplus <- rep(0,nsample)
    Zminus<- rep(0,nsample)
    
    mu1 <-  as.matrix(rbind(Covarmainphi)) %*% t(t(betaphi))
    mu2 <-  as.matrix(rbind(Covarmainmu)) %*% t(t(betamu))
    
    if (is.vector(Covarplus)) {Covarplusall <- c(Covarplus, Covarvalplus)} else {
      names(Covarvalplus) <- names(Covarplus) 
      Covarplusall <- rbind(Covarplus, Covarvalplus)
    }
    
    if (is.vector(Covarminus)) {Covarminusall <- c(Covarminus, Covarvalminus)} else {
      names(Covarvalminus) <- names(Covarminus) 
      Covarminusall <- rbind(Covarminus, Covarvalminus)
    }
    
    muZplus <- as.matrix(Covarplusall) %*% t(t(alphaplus)) 
    muZminus <- as.matrix(Covarminusall) %*% t(t(alphaminus)) 
    U1 <- rep(0,nmain)
    U2 <- rep(0,nmain)
    
    Y <- Ystar - Zplus[1:nmain] + Zminus[1:nmain]
    Y <- c(Y,Yval)
    
    beta1record <- betaphi
    beta2record <- betamu
    alphaplusrecord <- alphaplus
    alphaminusrecord <- alphaminus
    
    for (repk in 1:nmcmc){
      
      if (repk %% 500 ==0){
        cat("MCMC interations:" , repk,"\n")
        }
      
      ### Step 1: Measurement Error: Positive Process
      
      for (t in 1:nmain){
        NewGen <-  ZI_GenerateBigJoint(Ystar[t], Uibound[1], Uibound[2],  mu1[t], mu2[t], muZplus[t],  muZminus[t])
        Zplus[t] <- NewGen[1]
        Zminus[t] <- NewGen[2]
        Y[t] <- NewGen[3]
        U1[t] <- NewGen[4]
        U2[t] <- NewGen[5]
      }
      
      for (t in 1:nval){
        NewGen2 <- ZI_GenerateZpZmJoint(Ystarval[t], Yval[t], muZplus[t+nmain], muZminus[t+nmain])
        Zminus[t+nmain] <- NewGen2[1]
        Zplus[t+nmain] <- NewGen2[2]
      }
      
      alphaplussafe <- alphaplus
      
      alphaplus <- ZI_GenerateBetaMetro(alphaplus, Zplus, as.matrix(Covarplusall), propsigmaplus,  priorgamma)
      
      muZplus <- as.matrix(Covarplusall) * alphaplus
      
      
      # ### Step 2: Measurement Error: Negative Process
      
      alphaminussafe <- alphaminus
      
      alphaminus <- ZI_GenerateAlphaMNMetro(alphaminus, Y, Zminus, as.matrix(Covarminusall),  propsigma= propsigmaminus,  priormu,  priorSigmas=priorSigma)
      
      muZminus <- as.matrix(Covarminusall) * alphaminus
      
      # ### Step 3: Main Model
      
      ### Step 3.1 Update betaphi and beta 2
      
      betaphi <- ZI_GenerateBetaMetro(Par=betaphi, U=U1, Covar=as.matrix(Covarmainphi), propsigma=propsigmaphi, priorgamma)
      mu1 <-  as.matrix(Covarmainphi) %*% t(t(betaphi))
      
      betamu <- ZI_GenerateBetaMetro(Par=betamu, U=U2, Covar=as.matrix(Covarmainmu), propsigma=propsigmamu, priorgamma)
      
      mu2 <- as.matrix(Covarmainmu) %*% t(t(betamu))
      
      
      beta1record <- rbind(beta1record,betaphi)
      beta2record <- rbind(beta2record,betamu)
      if (is.vector(Covarplus)) {alphaplusrecord <- c(alphaplusrecord,alphaplus)} else {alphaplusrecord <- rbind(alphaplusrecord,alphaplus)}
      if (is.vector(Covarminus)) {alphaminusrecord <- c(alphaminusrecord,alphaminus)} else {alphaminusrecord <- rbind(alphaminusrecord,alphaminus)}
    }
    
    BayesResults <- list(betamu_trace = beta1record, betaphi_trace = beta2record,  
                         alphaplus_trace = alphaplusrecord, alphaminus_trace = alphaminusrecord)
    
    class(BayesResults) <- "ZIPBayes"
    
    
    return(BayesResults)

}

