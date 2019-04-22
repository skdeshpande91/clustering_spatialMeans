path <- "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/"

library(Rcpp)
library(RcppArmadillo)



source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

sourceCpp("src/ensm_cluster_mean.cpp")
sourceCpp("src/map_cluster.cpp")
sourceCpp("src/local_search.cpp")


load(paste0(path, "data/small_example_1.RData"))


map0_small <- map_partition(Y = y_1_small, A_block = A_block_small, gamma = gamma_0_small, a1 = 1/10, a2 = 10)
particle0_small <- ensm_cluster_mean(Y = y_1_small, A_block = A_block_small, L = 25, gamma = gamma_0_small, a1 = 1/10, a2 = 10, lambda = 100)
local0_small <- local_search(Y = y_1_small, A_block = A_block_small, gamma = map0_small$particles[[1]], a1 = 1/10, a2 = 10)

map2_large <- map_partition(Y = y_1_large, A_block = A_block_large, gamma = gamma_2_large, a1 = 1/10, a2 = 10)
particle2_large <- ensm_cluster_mean(Y = y_1_large, A_block = A_block_large, L = 5, gamma = gamma_2_large, a1 = 1/10, a2 = 10, lambda = 400)
local2_large <- local_search(Y = y_1_large, A_block = A_block_large, gamma = map2_large$particles[[1]], a1 = 1/10, a2 = 10)


test_map <- map_partition(Y = y_1_small, A_block = A_block_small, gamma = gamma_2_small, a1 = 1/10, a2 = 10)
test_particle <- ensm_cluster_mean(Y = y_1_small, A_block = A_block_small, L = 25, gamma_init = gamma_2_small, a1 = 1/10, a2 = 10, lambda = 100)
gamma_map <- test_map$map[[1]]



# why does the MAP procedure stop when the particle procedure does not.

test_particle <- ensm_cluster_mean(Y = y_1_small, A_block = A_block_small, L = 5, gamma_init = gamma_2_small, a1 = 1/10, a2 = 10, lambda = 1)

test_particle <- ensm_cluster_mean(Y = y_1_small, A_block = A_block_small, L = 5, gamma_init = gamma_1_small,
                                   a1 = 1/10, a2 = 10, lambda = 100)

test_map_large <- map_partition(Y = y_1_large, A_block= A_block_large, gamma = gamma_2_large, a1 = 1/10, a2 = 10)
test_particle_large <- ensm_cluster_mean(Y = y_1_large, A_block = A_block_large, L = 5, gamma_init = gamma_2_large, 
                                         a1 = 1/10, a2 = 10, lambda = 100)
