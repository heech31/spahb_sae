model.lik <- function(theta,xx,beta.post,Om.post){

  cc <-  as.vector( theta - xx%*%beta.post )
  qq <- sum( cc*rowSums( Om.post*cc ) )
  val <- as.numeric( 0.5*determinant(Om.post)$modulus - qq/2 )
  return(val)

}
