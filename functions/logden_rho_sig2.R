# Marginal log density of rho and sigma2

logden <- function(para, W, type=NULL){
  # para[1] = sigma2 : Model error variance parameter
  # para[2] = rho    : Spatial autocorrelation parameter
  # W      : Adjacency matrix
  # type   : type of W (FH,CAR,SAR,IAR,SRM)
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  sigma2 <- para[1]
  
  rho <- para[2]
  
  V   <-  D + solve(Omega(sigma2,rho,W, type=type))
  
  iVX <- solve(V,X)
  
  iVy <- as.vector( solve(V,Y) )
  
  XiVX <- crossprod( X, iVX )
  
  Q      <- quad.form( as.matrix( V - X %*% solve(XiVX,t(X))  ),  iVy )
  
  logdet <- determinant(V)$modulus + determinant(XiVX)$modulus
  
  logden <- as.numeric( -0.5*(logdet + Q) )
  
  return( logden )
}









logden.contour <- function(sigma2, rho, W, type=NULL){
  # para[1] = sigma2 : Model error variance parameter
  # para[2] = rho    : Spatial autocorrelation parameter
  # W      : Adjacency matrix
  # type   : type of W (FH,CAR,SAR,IAR,SRM)
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  V   <-  D + solve(Omega(sigma2,rho,W, type=type))
  
  iVX <- solve(V,X)
  
  iVy <- as.vector( solve(V,Y) )
  
  XiVX <- crossprod( X, iVX )
  
  Q      <- quad.form( as.matrix( V - X %*% solve(XiVX,t(X))  ),  iVy )
  
  logdet <- determinant(V)$modulus + determinant(XiVX)$modulus
  
  logden <- as.numeric( -0.5*(logdet + Q) )
  
  return( logden )
}













