library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

sourceCpp(file = "src/local_search.cpp")

load("data/partitions_small.RData")
load("data/small_example_1.RData")
load("data/scaled_alpha_small.RData")


small_0 <- local_search(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, eta = 1)
small_1 <- local_search(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 1)
small_2 <- local_search(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_2_small, a1 = 1/10, a2 = 10, eta = 1)
small_3 <- local_search(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_3_small, a1 = 1/10, a2 = 10, eta = 1)
small_4 <- local_search(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_4_small, a1 = 1/10, a2 = 10, eta = 1)
save(small_0, small_1, small_2, small_3, small_4, file = "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/local_search_small.RData")

# Larger grid
load("data/partitions_large.RData")
load("data/large_example_1.RData")
load("data/scaled_alpha_large.RData")
large_0 <- local_search(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large, a1 = 1/10, a2 = 10, eta = 1)
large_1 <- local_search(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, eta = 1)
large_2 <- local_search(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_2_large, a1 = 1/10, a2 = 10, eta = 1)
large_3 <- local_search(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_3_large, a1 = 1/10, a2 = 10, eta = 1)
large_4 <- local_search(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_4_large, a1 = 1/10, a2 = 10, eta = 1)
save(large_0, large_1, large_2, large_3, large_4, file = "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/local_search_large.RData")

