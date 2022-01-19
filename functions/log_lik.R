library(mvtnorm)

log_lik <- function(Y,theta.post,sDi){
  
  val <- dnorm(Y,theta.post,sDi, log=TRUE)
  
  return(val)
  
}

