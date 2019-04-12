library(Rcpp)
library(RcppArmadillo)
load("data/partitions_small.RData")
load("data/small_example_1.RData")
load("data/scaled_alpha_small.RData")
source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

sourceCpp(file = "src/local_search.cpp")

test_local <- local_search(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_0_small, 
                           a1 = 1/10, a2 = 10, eta = 1)
