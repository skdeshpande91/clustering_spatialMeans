# Test the new code for finding spectral splits

library(Rcpp)
library(RcppArmadillo)


source("scripts/partition_functions.R")
source("scripts/plot_partition.R")
load("data/partitions_large.RData")
load("data/large_example_1.RData")

sourceCpp("src/test_splits.cpp")


test <- test_moves(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_1_large, a1 = 1/10, a2 = 10)

test <- test_moves(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_0_large, a1 = 1/10, a2 = 10)


islands <- test_island(ybar_1_large[,1], T = 10, A_block_large, gamma_0_large, a1 = 1/10, a2 = 10, island_frac = 0.01)

borders <- test_border(ybar_1_large[,1], T = 10, A_block_large, gamma_0_large, a1 = 1/10, a2 = 10, rho = 0.99)


test_tail <- test_tail_split(ybar_1_large[,1], T = 10, A_block_large, gamma_1_large, a1 = 1/10, a2 = 10, rho = 0.99)



test <- test_spec_split(ybar_1_large[,1], T = 10, A_block_large, gamma_1_large, a1 = 1/10, a2 = 10, rho = 0.99, reps = 1000)


################
# Testing for main function 17:34 13 March 2019
sourceCpp("src/ensm_cluster_mean2.cpp")
test <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10, gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1,lambda = 1)

