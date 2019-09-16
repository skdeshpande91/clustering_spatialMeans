# Test the new ensm_cluster_mean
setwd("~/Documents/clustering_spatialMeans/")
library(Rcpp)
library(RcppArmadillo)



source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/partition_functions.R")
#source("~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/scripts/plot_partition.R")
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

set.seed(seed_seq[sim_number] * 100 + 10*(batch - 1) + r)
Y <- matrix(nrow = N, ncol = T)

for(i in 1:N) Y[i,] <- rnorm(T, mean = alpha[i, r + (batch-1)*batch_size], sd = 1)

sourceCpp("src/ensm_cluster_mean.cpp")
lam_1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, rho = 0.9)
