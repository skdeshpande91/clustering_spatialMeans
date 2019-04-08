library(Rcpp)
library(RcppArmadillo)
load("data/partitions_small.RData")
load("data/small_example_1.RData")
load("data/scaled_alpha_small.RData")

load("data/partitions_large.RData")
load("data/large_example_1.RData")
load("data/scaled_alpha_large.RData")

source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

sourceCpp("src/map_cluster.cpp")
sourceCpp("src/ensm_cluster_mean.cpp")

test0_small <- map_partition(ybar_1_small[,1], T = 10, A_block_small, gamma_init = gamma_0_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle0_small <- ensm_cluster_mean(ybar_1_small[,1], T = 10, L = 10, A_block_small, gamma_init = gamma_0_small, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)

test1_small <- map_partition(ybar_1_small[,1], T = 10, A_block_small, gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle1_small <- ensm_cluster_mean(ybar_1_small[,1], T = 10, L = 10, A_block_small, gamma_init = gamma_1_small, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)

test2_small <- map_partition(ybar_1_small[,1], T = 10, A_block_small, gamma_init = gamma_2_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle2_small <- ensm_cluster_mean(ybar_1_small[,1], T = 10, L = 10, A_block_small, gamma_init = gamma_2_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)

test3_small <- map_partition(ybar_1_small[,1], T = 10, A_block_small, gamma_init = gamma_3_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle3_small <- ensm_cluster_mean(ybar_1_small[,1], T = 10, L = 10, A_block_small, gamma_init = gamma_3_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)


test4_small <- map_partition(ybar_1_small[,1], T = 10, A_block_small, gamma_init = gamma_4_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle4_small <- ensm_cluster_mean(ybar_1_small[,1], T = 10, L = 10, A_block_small, gamma_init = gamma_4_small, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)

save(test0_small, particle0_small, test1_small, particle1_small, test2_small, particle2_small, test3_small, particle3_small, test4_small, particle4_small, 
     file = "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/test_map_small_results.RData")


# Now try large

test0_large <- map_partition(ybar_1_large[,1], T = 10, A_block_large, gamma_init = gamma_0_large, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle0_large <- ensm_cluster_mean(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_init = gamma_0_large, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)

test1_large <- map_partition(ybar_1_large[,1], T = 10, A_block_large, gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
# test1_large gets stuck at initial point but particle1_large does not
particle1_large <- ensm_cluster_mean(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_init = gamma_1_large, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)

test2_large <- map_partition(ybar_1_large[,1], T = 10, A_block_large, gamma_init = gamma_2_large, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle2_large <- ensm_cluster_mean(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_init = gamma_2_large, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)

test3_large <- map_partition(ybar_1_large[,1], T = 10, A_block_large, gamma_init = gamma_3_large, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle3_large <- ensm_cluster_mean(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_init = gamma_3_large, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)

test4_large <- map_partition(ybar_1_large[,1], T = 10, A_block_large, gamma_init = gamma_4_large, a1 = 1/10, a2 = 10, eta = 1, max_iter = 20)
particle4_large <- ensm_cluster_mean(ybar_1_large[,1], T = 10, L = 10, A_block_large, gamma_init = gamma_4_large, a1 = 1/10, a2  = 10 ,eta = 1, max_iter = 20)



save(test0_large, particle0_large, test1_large, particle1_large, test2_large, particle2_large, test3_large, particle3_large, test4_large, particle4_large,
     file = "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/test_map_large_results.RData")
