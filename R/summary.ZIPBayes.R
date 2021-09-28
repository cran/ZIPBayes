summary.ZIPBayes <-function (object, burnin = 1, thinperiod = 1, confidence.level = 0.95, ...) 
{ 
  
  confidence.margin <-  (1-confidence.level)/2
  
  summarized.res <- lapply( object, FUN= function(resultdata){
    results <- apply(resultdata, MARGIN =2, FUN = function(x){
      xpure <- purifyseq(x, burnin, thinperiod)
      mean <- mean(xpure,na.rm=T)
      sd <- sd(xpure,na.rm=T)
      median <- median(xpure,na.rm=T)
      CrI <- quantile(xpure, c(confidence.margin, 1 - confidence.margin),na.rm=T)
      HDR <- getHDP(xpure,confidence.margin*2)
      return(c(mean,median,sd,CrI,HDR))
    })
    rownames(results) <- c("mean", "sd", "median", paste0("CI:",confidence.margin,"%"), paste0("CI:", 1 - confidence.margin,"%"),
             "HDR_LB","HDR_UB")
    return(results)
  })
  

  class(summarized.res) <- "summary.ZIPBayes"
  
  summarized.res
  
}
