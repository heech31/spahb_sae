rho.bound <- function(W, type=NULL){
  
  if(is.null(type)){print("Error: Spatial structure needs to be specified"); break}
  
  if( type == "fh"){# FH model
    
    lk <- 0; uk <- 1;
    
  }else if( type == "car"){# CAR
    
    lambda <- eigen(W)$values
    
    lk <- 1/min(lambda); uk <- 1/max(lambda);
    
  }else if ( type =="sar"){
    
    lk <- -1; uk <- 1;
      
  }else if ( type == "iar"){
    
    lk <- -1; uk <- 1;
      
  }else if (type == "srm"){
    
    lk <- 0; uk <-  1;
    
  }
 return(c(lk,uk)) 
}