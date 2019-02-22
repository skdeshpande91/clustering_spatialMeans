library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/plot_partition.R")


load("data/partitions_large.RData")
load("data/large_example_1.RData")

sourceCpp("src/ensm_cluster_mean2.cpp")

test_1_L10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L = 10, gamma_init = gamma_1_large, a1 = 1/20, a2 = 20, max_iter = 20, eta = 1, split_frac = 0.2)


test_0_L10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L= 10, gamma_init = gamma_0_large, a1 = 1/20, a2 = 20, max_iter = 20, eta = 1)
plot_partition_grid(test_1_L10$particles[[4]], A_block_large, test_1_L10$alpha_hat_particle[,4])


test_2_L10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L = 10, gamma_init = gamma_2_large, a1 = 1/20, a2 = 20, max_iter = 20, eta = 1)


test_1_L10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L = 10, gamma_init = gamma_1_large, a1 = 1/20, a2 = 20, max_iter = 20, eta = 1, split_frac = 0.2)


test_0_L20 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L= 20, gamma_init = gamma_0_large, a1 = 1/20, a2 = 20, max_iter = 20, eta = 5)





test_1_L10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L = 10, gamma_init = gamma_1_large, a1 = 1/20, a2 = 20)
test_1_L10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block_large, L = 20, gamma_init = gamma_1_large, a1 = 1/20, a2 = 20)
