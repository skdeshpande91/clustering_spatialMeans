# Testing large example
library(Rcpp)
library(RcppArmadillo)
load("data/partitions_large.RData")
load("data/large_example_1.RData")

sourceCpp("src/ensm_cluster_mean.cpp")
test_0 <- ensm_cluster_mean(ybar_1_large[,1], T = 10, A_block_large, 
                            L= 10, gamma_init = gamma_0_large, a1 = 1/20, a2 = 20, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)

test_1 <- ensm_cluster_mean(ybar_1_large[,1], T = 10, A_block_large, 
                            L= 10, gamma_init = gamma_1_large, a1 = 1/20, a2 = 20, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)
