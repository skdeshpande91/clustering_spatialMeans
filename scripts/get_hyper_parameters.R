get_hyper_parameters <- function(Y, K, rho, hp_id = c(0, 1, 2, 3)){
  
  
  N <- nrow(Y)
  T <- ncol(Y)
  alpha_mle <- rowMeans(Y)
  sigma2_mle <-  (T-1)/T* apply(Y, MAR = 1, FUN = var)
  
  mean_sigma2 <- mean(sigma2_mle)
  var_sigma2 <- var(sigma2_mle)
  
  nu_sigma <- 2*((mean_sigma2 * mean_sigma2)/var_sigma2 + 2)
  lambda_sigma <- mean_sigma2 * (nu_sigma-2)/nu_sigma
  sigma2_est <- mean(sigma2_mle)
  
  if(hp_id == 0){
    hyper_params <- list(a1 = 1/T, a2 = T, nu_sigma = 3, lambda_sigma = 1)
  } else if(hp_id == 1){
    hyper_params <- list(a1 = 1/T, a2 = T, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma)
  } else if(hp_id == 2){
    a2 <- (max(abs(alpha_mle)))^2/(4 * sigma2_est)
    sigma_cl <- (max(alpha_mle) - min(alpha_mle))/(2 * (K+1))
    a1 <- (1 - rho) * sigma_cl * sigma_cl/(sigma2_est)
    hyper_params <- list(a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma)
  } else if(hp_id == 3){
    sigma_cl <- (max(alpha_mle) - min(alpha_mle))/(2 * (K+1))
    a1 <- (1 - rho) * sigma_cl * sigma_cl/(sigma2_est)
    a2 <- max(abs(alpha_mle))^2/(4 * sigma2_est) - a1/(1 - rho)
    if(a2 < 0) hyper_params <- NULL
    else hyper_params <- list(a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma)
  }

  return(hyper_params)
  
  
}
