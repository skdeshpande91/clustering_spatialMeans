# Quick script to compare different initializations
setwd("~/Documents/clustering_spatialMeans/")
library(Rcpp)
library(RcppArmadillo)




source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/partition_functions.R")
source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/plot_partition.R")
source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/get_hyper_parameters.R")


load("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/data/alphas.RData")
load("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/data/partitions.RData")



alpha <- alpha_1
N <- 400
T <- 12
Y <- matrix(nrow = N, ncol = T)
sim_number <- 1
r <- 1
batch <- 1
batch_size <- 5
seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)
hyper_param_id <- 3


set.seed(seed_seq[sim_number] * 100 + 10*(batch - 1) + r)
Y <- matrix(nrow = N, ncol = T)
for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i, r + (batch-1)*batch_size], sd = 1)
hyper_params <- get_hyper_parameters(Y, floor(log(N)), rho = 0.9, hp_id = 3)

a1 <- hyper_params[["a1"]]
a2 <- hyper_params[["a2"]]
nu_sigma <- hyper_params[["nu_sigma"]]
lambda_sigma <- hyper_params[["lambda_sigma"]]

sourceCpp("src/ensm_cluster_mean.cpp")

test_init0_lam1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 0, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)
test_init1_lam1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 1, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)
test_init2_lam1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 2, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)

save(test_init0_lam1, test_init1_lam1, test_init2_lam1, file = "alpha1_init_test_lam1.RData")



test_init0_lam100 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 0, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)
test_init1_lam100 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 1, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)
test_init2_lam100<- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 2, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)

save(test_init0_lam100, test_init1_lam100, test_init2_lam100, file = "alpha1_init_test_lam100.RData")


## Make a plot showing off the different initializations

png("images/alpha1_init0_lam100.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,4))
for(i in 1:10){
  plot_partition_grid(test_init0_lam100$particles[[i]], W, values = test_init0_lam100$alpha_hat_particle[,i], max_value = max_value,
                      title = paste0("log_post = ", round(test_init0_lam100$log_post[i], digits = 4)))
}
dev.off()

png("images/alpha1_init1_lam100.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,4))
for(i in 1:10){
  plot_partition_grid(test_init1_lam100$particles[[i]], W, values = test_init1_lam100$alpha_hat_particle[,i], max_value = max_value,
                      title = paste0("log_post = ", round(test_init1_lam100$log_post[i], digits = 4)))
}
dev.off()

png("images/alpha1_init2_lam100.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,4))
for(i in 1:10){
  plot_partition_grid(test_init2_lam100$particles[[i]], W, values = test_init2_lam100$alpha_hat_particle[,i], max_value = max_value,
                      title = paste0("log_post = ", round(test_init2_lam100$log_post[i], digits = 4)))
}
dev.off()




### Now let's do something less well-separated
sim_number <- 4
alpha <- get(paste0("alpha_", sim_number))

r <- 1
batch <- 1
batch_size <- 5
seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)
hyper_param_id <- 3

set.seed(seed_seq[sim_number] * 100 + 10*(batch - 1) + r)
Y <- matrix(nrow = N, ncol = T)
for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i, r + (batch-1)*batch_size], sd = 1)
hyper_params <- get_hyper_parameters(Y, floor(log(N)), rho = 0.9, hp_id = 3)

test_init0_lam1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 0, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)
test_init1_lam1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 1, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)
test_init2_lam1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 2, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 1)

save(test_init0_lam1, test_init1_lam1, test_init2_lam1, file = "alpha4_init_test_lam1.RData")



test_init0_lam100 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 0, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)
test_init1_lam100 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 1, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)
test_init2_lam100<- ensm_cluster_mean(Y = Y, A_block = W, L = 10, init_id = 2, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9, lambda = 100)

save(test_init0_lam100, test_init1_lam100, test_init2_lam100, file = "alpha4_init_test_lam1.RData")

png("images/alpha4_init0_lam100.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,4))
for(i in 1:10){
  plot_partition_grid(test_init0_lam100$particles[[i]], W, values = test_init0_lam100$alpha_hat_particle[,i], max_value = max_value,
                      title = paste0("log_post = ", round(test_init0_lam100$log_post[i], digits = 4)))
}
dev.off()

png("images/alpha4_init1_lam100.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,4))
for(i in 1:10){
  plot_partition_grid(test_init1_lam100$particles[[i]], W, values = test_init1_lam100$alpha_hat_particle[,i], max_value = max_value,
                      title = paste0("log_post = ", round(test_init1_lam100$log_post[i], digits = 4)))
}
dev.off()

png("images/alpha4_init2_lam100.png", width = 6, height = 4.5, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,4))
for(i in 1:10){
  plot_partition_grid(test_init2_lam100$particles[[i]], W, values = test_init2_lam100$alpha_hat_particle[,i], max_value = max_value,
                      title = paste0("log_post = ", round(test_init2_lam100$log_post[i], digits = 4)))
}
dev.off()







