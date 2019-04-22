# Testing small

# Run the particle search for small grid examples


library(Rcpp)
library(RcppArmadillo)
load("data/partitions_small.RData")
load("data/small_example_1.RData")
load("data/scaled_alpha_small.RData")
source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

sourceCpp("src/ensm_cluster_mean.cpp")


test0 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 5, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20, lambda = 100)
test1 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 5, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20, lambda = 100)
test2 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 5, gamma_init = gamma_2_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20, lambda = 100)

plot_particle_set_new(test1$particle_trajectory[[5]], A_block_small)


sourceCpp("src/ensm_cluster_mean2.cpp")


test1 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 5, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20, lambda = 100)
test2 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20, lambda = 100)


test1 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
test2 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)


test1 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 5, max_iter = 20)
test2 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 5, max_iter = 20)



test_1 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.2, max_iter = 20)

test2 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)


test_3 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, L = 50, gamma_init = gamma_3_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)

test_0 <- ensm_cluster_mean(ybar_1_small[,1], T = 10, A_block_small, 
                            L= 50, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, rho = 0.99, eta = 1, lambda = 1, split_frac = 0.1, max_iter = 20)
