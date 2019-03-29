# testing progressive merges
library(Rcpp)
library(RcppArmadillo)
load("data/partitions_small.RData")
load("data/small_example_1.RData")

sourceCpp("src/test_prog_merge.cpp")

test <- test_prog_merge(ybar_1_small[,1], T = 10, A_block = A_block_small, gamma_init = gamma_0_small, L = 10)
