# Sampling theta from its conditional posterior distribution
theta_sample <-  function(beta,sigma2,rho,W,type=NULL,X){
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  Om   <- as.matrix( Omega(sigma2=sigma2, rho=rho, W=W, type=type) )
  
  Dinv <- diag(1/diag(D))
  
  theta_Sigma <- solve(Dinv + Om )
  
  theta_mu    <-  Y - theta_Sigma%*%Om%*%(Y-X%*%beta)
  
  theta <- rmvnorm(1, mean=theta_mu, sigma=theta_Sigma)
  
  return(theta)
  
}
