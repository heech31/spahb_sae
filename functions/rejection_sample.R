rejection.sampling <- function(Y, X, D, W, type=NULL, S, verbose=FALSE){
  
  
  # Y : m x 1 with direct estimates
  # X : m x p covariates (auxiliary variables) in columns. Assume to be standardized
  # D : m x m diagonal matrix with sampling variances
  # type : Spatial structure (FH, CAR, SAR, IAR, SRM)
  # S    : Posterior sample size
  
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  m <- nrow(X) # Number of small areas
  p <- ncol(X) # Number of covariates including intercept

  lkuk <- rho.bound(W,type) # lower (lk) and upper (uk) bounds of rho
  lk <- lkuk[1]+1e-4; uk <- lkuk[2]-1e-4
  
  # find the mode of the (unnormalized) log density of (sigma2, rho)
  if(type!="fh"){
    
    optimized <- optim(c(1,0.5),function(x) logden(x,W,type=type),
                       lower=c(1e-8,lk),upper=c(1e+8,uk), method = "L-BFGS-B", 
                       control=list(fnscale=-1))
    # Mode 
    optim.mean  <- optimized$par 
    # -Hessian inverse
    optim.sigma <- 2*(-solve( hessian(function(x) logden(x,W,type=type), optimized$par) ))
  
  }else{
    
    optimized <- optimize(function(x) logden(c(x,0),W,type=type), lower=0, upper=100, maximum=TRUE )
    
    # Mode 
    optim.mean  <- optimized$maximum
    # -Hessian inverse
    optim.sigma <- 2*(-solve( hessian(function(x) logden(c(x,0),W,type=type), optim.mean) ))
  
  }
  
  
  if( any(eigen(optim.sigma)$values<0) ){
    
    optim.sigma <- diag(abs(diag(optim.sigma)) ) + diag(2)*0.1
  }
  
  if( any(diag(optim.sigma)<1e-2) ){
    
    diag(optim.sigma)[ diag(optim.sigma)<1e-3 ] <- 0.1
    
  }
  
  # log(f/g) to find M
  tmpfun <- function(x){ 
    fun.val <- logfg(x, W=W, lk=lk, uk=uk, type=type, 
                     optim.mean = optim.mean, optim.sigma = optim.sigma )
    return(fun.val)
           }

  
  # log(M) for log( f/(Mg) ) = log(f) - log(M) - log(g)
  if( type != "fh"){
    
    logM <- optim(c(1,0.5), tmpfun, method = "L-BFGS-B", 
                  lower=c(1e-8,lk),upper=c(1e+8,uk), 
                  control=list(fnscale=-1))$value
  }else{
    
    logM <- optimize(f = tmpfun, interval = c(1e-8,1e+8), maximum=TRUE)$objective
  
  }
  
  
  # Number of proposals  ntry/S = acceptance rate
  ntry <- 0

  # Posterior sample of sigma2 and rho
  post_s2_rho <- matrix(NA,2,S)
  # Posterior sample of beta
  post_beta  <- matrix(NA,p,S)
  # Posterior sample of theta
  post_theta <- matrix(NA,m,S)

  bt <- Sys.time()
  break.flag <- FALSE
  ################################################################################
  # Rejection sampling
  if( type != "fh"){
    
    for( ss in 1:S){
  	  lratio <- -Inf
  	  u     <- runif(1)
      while.count <- 0
  	  while( lratio<= log(u) ){
  	  
  		  u     <- runif(1)
  		  x     <- as.vector( tmvtnorm::rtmvt(1, mean=optim.mean, sigma = optim.sigma, 
  		                            df=3, lower=c(1e-8,lk), upper=c(1e+8,uk) ) )
  
  		  logd_target <- logden(para=x, W=W, type=type)
  		  logd_prop   <- tmvtnorm::dtmvt(x, mean=optim.mean, sigma = optim.sigma, 
  		                       df=3, lower=c(1e-8,lk), upper=c(1e+8,uk), log=TRUE)
  		  lratio      <- logd_target - logd_prop - logM
  		  ntry        <- ntry + 1
  		  #print(lratio<= logM+log(u))
  		
  		  while.count <- while.count + 1
  		  if(while.count%%1000==0){
  		    print("Acceptance rate is too low")
  		    break.flag <- TRUE
  		    break;
  		  }
  		  if(break.flag){break;}
  	  }
      

  	  post_s2_rho[,ss] <- x
  	  if(verbose==TRUE){
  	    if( ss%%( round(S*0.05) ) == 0 ) print(ss/ntry)
  	    }
  	  
  	  if(break.flag){break;}
  	  
  	  }# lratio-logM
    
  }else{
    
    for( ss in 1:S){
      lratio <- -Inf
      u     <- runif(1)
      while.count <- 0
      while( lratio<= log(u) ){
        
        u     <- runif(1)
        x     <- as.vector( TruncatedNormal::rtmvt(1, mu=optim.mean, sigma = optim.sigma, 
                                  df=3, lb=c(0), ub=c(Inf) ) )
        
        logd_target <- logden(para=c(x,0), W=W, type=type)
        logd_prop   <- TruncatedNormal::dtmvt(x, mu=optim.mean, sigma = optim.sigma, 
                               df=3, lb=c(0), ub=c(Inf), log=TRUE)
        lratio      <- logd_target - logd_prop - logM
        ntry        <- ntry + 1
        
        while.count <- while.count + 1
        if(while.count%%1000==0){
          print("Acceptance rate is too low")
          print(while.count)
        }
        
      }
      post_s2_rho[1,ss] <- x
      if(verbose==TRUE){
        if( ss%%( round(S*0.05) ) == 0 ) print(ss/ntry)
      }
    }# lratio-logM
    
    
  }
    
  #
  ################################################################################
  
  for( ss in 1:S){
    if(break.flag){break;}
    post_beta[,ss] <- beta_sample(post_s2_rho[1,ss], 
                                  post_s2_rho[2,ss], W=W, type=type)
    post_theta[,ss] <- theta_sample(post_beta[,ss,drop=FALSE],post_s2_rho[1,ss], 
                                    post_s2_rho[2,ss], W=W, type=type,X)
    if(verbose==TRUE){
      if( ss%%( round(S*0.1) ) == 0 ) print(ss)
    }
  }
  ################################################################################

  et <- Sys.time()
  sampling.time <- et-bt
  acceptance.prob <- ss/ntry

  results <- list(
    post_s2_rho = post_s2_rho,
    post_beta = post_beta,
    post_theta = post_theta,
    sampling.time = sampling.time,
    acceptance.prob = acceptance.prob)
  
  
  return(results)

}
