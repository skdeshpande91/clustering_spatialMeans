library(Rcpp)
library(RcppArmadillo)


source("scripts/partition_functions.R")
source("scripts/plot_partition.R")
load("data/partitions_large.RData")
load("data/large_example_1.RData")


sourceCpp("src/test_km.cpp")

test <- test_km(ybar_1_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large, reps = 10)

