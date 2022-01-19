library(matrixcalc)
# Sampling beta from its conditional posterior distribution
beta_sample <- function(sigma2,rho,W,type=NULL){
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  V   <-  D + solve(Omega(sigma2,rho,W, type=type))
  
  VinvX <- solve(V,X)
  
  XVinvX <- crossprod( X, VinvX ) # X' V^-1 X
  
  XVinvX <- ( t(XVinvX) + XVinvX ) / 2
  
  Vinvy <- as.vector( solve(V,Y) ) # V^-1 Y

  Delta <- as.matrix( Matrix::solve( XVinvX ) )
  
  Delta <- (t(Delta) + Delta)/2
  
  gamma <- as.vector( Delta%*%t(X)%*%solve(V,Y) )
  
  beta_s <- rmvnorm(1,gamma,Delta)
  
  return(beta_s)
  
}  
