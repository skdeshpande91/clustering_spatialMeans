setwd("Documents/Research/clustering_spatialMeans/")

library(Rcpp)
library(RcppArmadillo)

source("scripts/partition_functions.R")
source("scripts/get_hyperparameters.R")



load("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/september2019/data/alphas.RData")
load("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/september2019/data/partitions.RData")

N <- 400
T <- 12

seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)

batch_size <- 5
batch <- 1
r <- 1

sim_number <- 1

alpha <- get(paste0("alpha_", sim_number))
set.seed(seed_seq[sim_number] * 100 + 10 * (batch - 1) + r)
Y <- matrix(nrow = N, ncol = T)
for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i, r + (batch-1)*batch_size], sd = 1)
hyper_params <- get_hyperparameters(Y, floor(log(N)), rho = 0.9)

a1 <- hyper_params[["a1"]]
a2 <- hyper_params[["a2"]]
nu_sigma <- hyper_params[["nu_sigma"]]
lambda_sigma <- hyper_params[["lambda_sigma"]]

sourceCpp("src/spectral_clustering.cpp")

test <- spectral_particle(Y, W, gamma_1, max_splits = 10, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)

