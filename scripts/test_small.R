# Testing small

# Run the particle search for small grid examples


library(Rcpp)
library(RcppArmadillo)
load("data/partitions_small.RData")
load("data/small_example_1.RData")

sourceCpp("src/ensm_cluster_mean.cpp")

test_1 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, 
                            L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)


test_3 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_3_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)

test_0 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, 
                            L= 50, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)
