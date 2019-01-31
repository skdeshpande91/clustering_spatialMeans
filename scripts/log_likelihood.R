# Log-likelihood for single sigma model
log_likelihood <- function(gamma, ybar, A_block, T = 10, a1 = 1/10, a2 = 10, rho = 0.99, a_sigma = 1, nu_sigma = 1)
{
  tmp <- partition_config(gamma, A_block)
  K <- tmp$K
  config <- tmp$config
  singletons <- tmp$singletons
  non_singletons <- tmp$non_singletons
  
  N <- length(ybar) 
  log_det_Omegay <- rep(NA, times = K)
  y_Omegay_y <- rep(NA, times = K)
  for(k in 1:K){
    if(config[k] == 1){
      tmp_Omegay = 1/(1/T + a2 + a1/(1- rho))
      log_det_Omegay[k] <- log(tmp_Omegay)
      y_Omegay_y[k] <- ybar[gamma[[k]]]^2 * tmp_Omegay
    } else{
      n_k <- config[k]
      y_k <- ybar[gamma[[k]]]
      A_block_k <- A_block[gamma[[k]], gamma[[k]]] # adjacency matrix of the cluster
      D <- diag(rowSums(A_block_k))
      A_star_k <- D - A_block_k # unweighted laplacian 
      Omega_alpha <- rho * A_star_k + (1 - rho) * diag(n_k)
      Sigma_alpha <- solve(Omega_alpha)
      Sigma_y <- a1 * Sigma_alpha + a2 + 1/T * diag(n_k)
      Omega_y <- solve(Sigma_y)
      #Sigma_y <- 1/a1 * Omega_alpha - (1/a1^2) * a2 * (1 - rho)^2/(1 + n_k/a1 * a2 * (1-rho))
      #diag(Sigma_y) <- diag(Sigma_y) + 1/T
      tmp_Omegay <- solve(Sigma_y)
      log_det_Omegay[k] <- determinant(tmp_Omegay, logarithm = TRUE)$modulus
      y_Omegay_y[k] <- t(ybar[gamma[[k]]]) %*% tmp_Omegay %*% ybar[gamma[[k]]]
    }
  }
  
  print("y_Omegay_y = ")
  print(y_Omegay_y)
  print("log_det = ")
  print(log_det_Omegay)
  ll <- 0.5 * sum(log_det_Omegay) - (a_sigma + N/2) * log( (nu_sigma + sum(y_Omegay_y))/2)
  return(ll)
}

alpha_postmean <- function(gamma, ybar, A_block, T = 10, a1 = 1/10, a2 = 10, rho = 0.99)
{
  tmp <- partition_config(gamma, A_block)
  K <- tmp$K
  config <- tmp$config
  singletons <- tmp$singletons
  non_singletons <- tmp$non_singletons
  
  N <- length(ybar) 
  alpha_hat <- rep(NA, times = N)
  for(k in 1:K){
    n_k <- config[k]
    if(n_k == 1){
      alpha_hat[gamma[[k]]] <- T * ybar[gamma[[k]]]/(T + 1/(a2 + a1/(1-rho)))
    } else{
      y_k <- ybar[gamma[[k]]]
      A_block_k <- A_block[gamma[[k]], gamma[[k]]] # adjacency matrix of the cluster
      D <- diag(rowSums(A_block_k))
      A_star_k <- D - A_block_k # unweighted laplacian 
      Omega_alpha <- rho * A_star_k + (1 - rho) * diag(n_k)
      V_inv <- 1/a1 * Omega_alpha - 1/a1^2 * (1-rho)^2/(1 + 1/a1 *a2 * n_k * (1-rho))
      diag(V_inv) <- diag(V_inv) + T
      V <- solve(V_inv)
      alpha_hat[gamma[[k]]] <- T * V %*% ybar[gamma[[k]]]
    }
  }
  return(alpha_hat)
}

