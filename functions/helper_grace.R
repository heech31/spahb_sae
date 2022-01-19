

loo_apply <-  function(Fitted){

  log_liks <- lapply(Fitted, function(x) extract_log_lik(x, merge_chains = FALSE))
  
  r_effs <- lapply(log_liks, function(x) relative_eff(exp(x)))
  
  loo.out <- Map(loo, log_liks, r_eff=r_effs, cores = getOption("mc.cores", 4))
  
  return(loo.out)

  }



