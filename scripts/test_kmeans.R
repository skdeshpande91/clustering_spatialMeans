setwd("~/Documents/Research/clustering_spatialMeans/")
library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/kmeans.cpp")

source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/partition_functions.R")
source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/plot_partition.R")
source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/get_hyper_parameters.R")


load("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/data/alphas.RData")
load("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/data/partitions.RData")


alpha <- alpha_1
N <- 400
T <- 12
sim_number <- 1
r <- 1
batch <- 1
batch_size <- 5
seed_seq <- c(129, 724, 603, 212, 1123, 2391, 815, 1947)
hyper_param_id <- 0

set.seed(seed_seq[sim_number] * 100 + 10*(batch - 1) + r)
Y <- matrix(nrow = N, ncol = T)
for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i, r + (batch-1)*batch_size], sd = 1)




hyper_params <- get_hyper_parameters(Y, floor(log(N)), rho = 0.9, hp_id = hyper_param_id)

a1 <- hyper_params[["a1"]]
a2 <- hyper_params[["a2"]]
nu_sigma <- hyper_params[["nu_sigma"]]
lambda_sigma <- hyper_params[["lambda_sigma"]]


km <- kmeans_particle(Y, W, gamma_1, max_split = 10, a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma, rho = 0.9)



test$init_particle
