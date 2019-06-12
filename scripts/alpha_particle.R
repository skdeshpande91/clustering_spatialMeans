# Functions used for the simulation setting
# This is not needed anymore

alpha_particle <- function(results, L = NULL){
  if(is.null(L)){
    L <- length(results$particles)
  }
  #weights <- results$pstar[1:L]
  weights <- results$pstar[1:L]/sum(results$pstar[1:L])
  tmp_alpha <- results$alpha_hat_particle
  alpha_hat <- rowSums(tmp_alpha %*% diag(weights))
  
  return(alpha_hat)
}

