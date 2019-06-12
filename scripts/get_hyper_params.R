#get_hyper_param <- function(Y, W, rho = 0.9, q = 0.05){
get_hyper_param <- function(Y, nu = 3,rho = 0.9, q = 0.05){
  #D <- diag(rowSums(W))
  #W_star <- D - W
  #Omega_alpha <- rho * W_star + (1 - rho) * diag(N)
  #Sigma_alpha <- solve(Omega_alpha)
  #avg_trace <- mean(diag(Sigma_alpha))
  
  T <- ncol(Y)
  alpha_hat <- rowMeans(Y)
  alpha_inf <- max(abs(alpha_hat))
  sigma2_alpha <- var(alpha_hat)
  sigma2_i <- apply(Y, FUN = var, MAR = 1) * (T-1)/T # don't want the unbiased estimates
  sigma2_hat <- mean(sigma2_i)
  
  lambda_sigma <- sigma2_hat * qchisq(0.1, df = nu)/nu
  if(alpha_inf^2/(sigma2_hat * qnorm(1 - q/2)^2) < sigma2_alpha/sigma2_hat){
    a2 <- alpha_inf^2/(sigma2_hat * qnorm(1 - q/2)^2)
    a1 <- (sigma2_alpha/sigma2_hat - a2) * (1 - rho)
    results <- c()
    results["lambda_sigma"] <- lambda_sigma
    results["a1"] <- a1
    results["a2"] <- a2
  } else{
    results <- NULL
  }
  return(results)
}
