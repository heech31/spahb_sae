generate.data <- function(m, sDi, sigma2, rho, WSAR, missing=NULL){

  A   <- solve(diag(m) - rho*(WSAR))

  vi  <- rnorm(m, 0, sqrt(sigma2) )

  x1star <- rnorm(m); x2star <- rnorm(m);

  x1  <- as.vector( A%*%x1star )

  x2	<- as.vector( A%*%x2star )

  th  <- 2.0*x1 + 1.0*x2  + vi

  y   <- th + rnorm(m,0,sDi)

  return( list(Y=y, X=cbind(1,x1,x2), theta=th ) )

}




