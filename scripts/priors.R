# Priors on gamma


# This is the original Ewens-Pitman prior
# leaving out the log(gamma(eta)) - log(gamma(eta + p)) from the mix here
# as it is constant, so long as eta is fixed.
log_pi_ep <- function(gamma, A_block, eta){
  tmp <- partition_config(gamma, A_block)
  K <- tmp$K
  config <- tmp$config
  
  #return(log(gamma(eta)) - log(gamma(eta + p)) + M * log(eta) + sum(log(gamma(config))))
  return(K*log(eta) + sum(lgamma(config)))
}

# How many clusters do we expect with ewens-pitman?

# Jensen-Liu prior
log_pi_jl <- function(gamma, A_block, eta){
  tmp <- partition_config(gamma, A_block)
  K <- tmp$K
  config <- tmp$config
  x <- sort(config) * log(seq(eta + 1, eta + K, length = K))
  return( (K - 1)*log(eta) + log(eta + K) + sum(x))
}

# Hierarchical Uniform prior
# we can use the hypergeo function / package to deal with the confluent hypergeometric function
# in the intrinsic prior
