rejection.sample.within.gibbs <- function(r.proposal, log.proposal, logM, log.target){
  
    r.proposal <- match.fun(r.proposal)
    log.proposal <- match.fun(log.proposal)
    log.target <- match.fun(log.target)
  
    lratio <- -Inf
    
    u     <- runif(1)
    
    while.count <- 0
    
    while( lratio<= log(u) ){
      
      u     <- runif(1)
      
      x     <- r.proposal()
      
      logd_target <- log.target(x)
      logd_prop   <- log.proposal(x)
      
      lratio      <- logd_target - logd_prop - logM #exp(lratio)
      ntry        <- ntry + 1
      #print(lratio<= logM+log(u))
      
      while.count <- while.count + 1
      if(while.count%%1000==0){
        print("Acceptance rate is too low")
        break.flag <- TRUE
        break;
      }
      
    }
    
    return(x)
}