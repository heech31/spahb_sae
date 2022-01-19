morans.I <- function(x,W){
  
  N <- length(x)
  
  x <-  x - mean(x)
  
  diag(W) <- 0
  
  tmp <- x%*%t(x)*W
  
  val <- ( sum(tmp)/sum( x^2 ) )*( N/sum( W ) )
  
  return(val)
  
}


# x <- d1989_cen[,"4pm"]/1000
morans.I.local <- function(x,W){
  
  N <- length(x)
  
  x <-  x - mean(x)
  
  diag(W) <- 0
  
  tmp <- x%*%t(x)*W
  
  val <- ( rowSums(tmp)/rowSums(W) )/ ( sum(x^2)/N )
  
  return(val)
  
}
