Omega <- function(sigma2,rho,W,type=NULL){
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  m  <- ncol(W)
  
  L <- rowSums(W)
  
  if( any(L==0) ){ L[L==0] <- 1}
  
  
  if( type == "fh"){# FH model
    
    Om <- diag(m)
    
    colnames(Om) <- rownames(Om) <- colnames(W)
    
    }else if( type == "car"){# CAR
      
      Om <- (diag(m) - rho*W)
      
    }else if ( type =="sar"){
      
      Wtild <- W*(1/L)
      
      Om_half <- (diag(m) - rho*Wtild)
      
      Om <- crossprod(Om_half)
      
      
    }else if ( type == "iar"){
      
      L <- diag(L)
      
      Om <- L - rho*W
      
    }else if (type == "srm"){
      
      L <- diag(L)
      
      R <- L - W
      
      Om <- rho*R + (1-rho)*diag(m)
      
    }
    
    return(Om/sigma2)
    
  }
  
