

logfg <- function(x, W, lk, uk, type=NULL, optim.mean, optim.sigma){
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  if( type != "fh" ){
    
    num  <- logden(x,W,type=type)
    
    denom<- tmvtnorm::dtmvt(x, mean=optim.mean, sigma = optim.sigma, df=3, lower=c(1e-8,lk), upper=c(1e+8,uk),log=TRUE) 
  
  }else{

    num  <- logden(c(x,0),W,type=type)
    
    denom<- TruncatedNormal::dtmvt(x, mu=optim.mean, sigma = optim.sigma, df=3, lb=c(1e-8), ub=c(1e+8),log=TRUE) 
    
  }
  return(num-denom)
}





