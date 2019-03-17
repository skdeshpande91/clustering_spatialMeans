library(Rcpp)
library(RcppArmadillo)


source("scripts/partition_functions.R")
source("scripts/plot_partition.R")
load("data/partitions_large.RData")
load("data/large_example_1.RData")


sourceCpp("src/test_km.cpp")

<<<<<<< HEAD
test <- test_km(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large, reps = 10)
=======

# Testing get_km_splits that has been updated with kmeans_repeat
# previously had observed running the basic test_km was returning a *different* set of splits
# than the one given by get_km_splits. The latter was always worse.

test_get_km <- test_km_particle(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10, gamma_1_large, a1 = 1/20, a2 = 20)

# Now try for a range of hyper-parameters values.



test_2 <- test_km(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large, reps = 500)



test_0 <- test_km_particle(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10, gamma_init = gamma_1_large, reps = 500)



test_1 <- test_kmpp(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large, reps = 500)


>>>>>>> 398318f1f16231f6f83f1f77d9e53c661affe695

