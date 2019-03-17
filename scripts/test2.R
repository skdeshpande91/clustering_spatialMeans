library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

load("data/partitions_large.RData")
load("data/large_example_1.RData")

load("data/partitions_small.RData")
load("data/small_example_1.RData")


sourceCpp("src/ensm_cluster_mean2.cpp")


#### Run small examples

small_L10_lam1 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 10, 
                              gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 1)
small_L10_lam5 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 10,
                                     gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 5)
small_L10_lam10 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 10,
                                     gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 5)

small_L50_lam1 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 50,
                                     gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 1)

small_L50_lam5 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 50,
                                     gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 5)

small_L50_lam10 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 50,
                                      gamma_init = gamma_1_small, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 10)
save(small_L10_lam1, small_L10_lam5, small_L10_lam10, file = "small_L10.RData")
save(small_L50_lam1, small_L50_lam5, small_L50_lam10, file = "small_L50.RData")


######## Run large examples
large_L10_lam1 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10,
                                     gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 1)
large_L10_lam5 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10,
                                     gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 5)
large_L10_lam10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10,
                                     gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 10)

save(large_L10_lam1, large_L10_lam5, large_L10_lam10, file = "large_L10.RData")

large_L25_lam1 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 25,
                                     gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 1)
large_L25_lam5 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 25,
                                     gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 5)
large_L25_lam10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 25,
                                      gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 10)
save(large_L25_lam1, large_L25_lam5, large_L25_lam10, file = "large_L25.RData")

# Run thre instances using the old version of best_split

test_1 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10, 
                           gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1,lambda = 1)


test_2 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 10, 
                             gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1,lambda = 1)

test_3 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10 ,A_block = A_block_large, L = 10,
                             gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 1)

test_25 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10 ,A_block = A_block_large, L = 25,
                              gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 1)


test_lambda10 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10 ,A_block = A_block_large, L = 25,
                              gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 10)
test_lambda5 <- ensm_cluster_mean2(ybar_1_large[,1], T = 10 ,A_block = A_block_large, L = 25,
                                   gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 5)

save(test_1, test_lambda5, test_lambda10, file = "test_results_mar15.RData")

test_1$pstar
test_1$counts
test_2$pstar
test_2$counts
test_3$pstar
test_3$counts







# A helpful plot
n_side <- 20
plot_partition_grid(gamma_0_large, A_block_large)
for(i in 1:n_side){
  for(j in 1:n_side){
    #rect(j-1, i-1, j, i, border = "lightgray", lty = 1, lwd = 0.5)
    text(j-0.5, i-0.5, labels = (j-1) + (i-1)*n_side, cex = 0.5) # stay consistent with C++ indexing
  }
}

# Run test if we start from the true partition
test_small_0 <- ensm_cluster_mean2(ybar_1_small[,1], T = 10, A_block = A_block_small, L = 10, gamma_init = gamma_0_small,
                                   a1 = 1/10, a2 = 10,nu_sigma = 3, lambda_sigma = 1)
