gibbs_rho_sig2 <- function(sig2, rho, lk, uk, burnin, nmc, verbose=FALSE){


  post_s2_rho <- matrix(NA,2,nmc)

  for( gg in 1:(burnin+nmc) ){
    
    log.sig2.given.rho <- Vectorize( function(x) logden(c(x,rho),W,type=type) )

    optim.mean <- optimize(log.sig2.given.rho, interval=c(1e-8,1e+8), maximum = TRUE)$maximum
    optim.var <- 1/(-hessian(log.sig2.given.rho, optim.mean))
    
    if(optim.var<1e-3){optim.var <- 0.01}

    log.proposal.sig2 <- Vectorize( function(x) 
      TruncatedNormal::dtmvt(x,mu=optim.mean,sigma=optim.var, df=4,lb=1e-8,ub=1e+8, log=TRUE) )

    sample.proposal.sig2 <- function(){
      TruncatedNormal::rtmvt(1,mu=optim.mean,sigma=optim.var,df=4,lb=1e-8,ub=1e+8 ) }

    objective <- function(x){ log.sig2.given.rho(x) - log.proposal.sig2(x)}

    logM <- optimize(objective,interval=c(1e-8,1e+8), maximum=TRUE)$objective
    
    
    sig2.new <- rejection.sample.within.gibbs(r.proposal = sample.proposal.sig2, 
                                         log.proposal = log.proposal.sig2,
                                         logM = logM,
                                         log.target = log.sig2.given.rho )

    sig2 <- sig2.new

###################################################################################################

    log.rho.given.sig2 <- Vectorize( function(x) logden(c(sig2,x),W,type=type) )

    optim.mean <- optimize(log.rho.given.sig2, interval=c(lk,uk), maximum = TRUE)$maximum
    optim.var <- 1/(-hessian(log.rho.given.sig2, optim.mean))
    
    if(optim.var<1e-3){optim.var <- 0.01}
    
    log.proposal.rho <- Vectorize( function(x) 
      TruncatedNormal::dtmvt(x,mu=optim.mean,sigma=optim.var,df=4,lb=lk,ub=uk, log=TRUE) )

    sample.proposal.rho <- function(){
      TruncatedNormal::rtmvt(1,mu=optim.mean,sigma=optim.var,df=4,lb=lk,ub=uk ) }
    
    objective <- function(x){ log.rho.given.sig2(x) - log.proposal.rho(x) }

    logM <- optimize(objective,interval=c(lk,uk), maximum=TRUE)$objective

    
        
    rho.new <- rejection.sample.within.gibbs(r.proposal = sample.proposal.rho, 
                                              log.proposal = log.proposal.rho,
                                              logM = logM,
                                              log.target = log.rho.given.sig2 )
    
    rho <- rho.new
    
    if(gg>burnin){
      post_s2_rho[,gg-burnin] <- c(sig2.new,rho.new)
    }
    
    if(verbose){
      if( gg%%round( (burnin+nmc)*0.1 ) == 0 ){
        print(gg)
        print( Sys.time() )
      }
    }
  }
  
  return(post_s2_rho)  

}
