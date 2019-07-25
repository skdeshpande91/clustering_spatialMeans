# Simulation for June 12
library(Rcpp)
library(RcppArmadillo)

dropbox_path <- "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/june5/"

load(paste0(dropbox_path, "data/partitions.RData"))
load(paste0(dropbox_path, "data/alphas.RData"))

source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

sourceCpp("~/Documents/clustering_spatialMeans/src/ensm_cluster_mean.cpp")
sourceCpp("~/Documents/clustering_spatialMeans/src/map_cluster.cpp")
sourceCpp("~/Documents/clustering_spatialMeans/src/partition_summary.cpp")

seed_seq <- c(129, 724, 603, 212)

for(sim_number in 1:4){
  print(paste("Starting sim number", sim_number, "at", Sys.time()))
  alpha <- get(paste0("alpha_", sim_number))
  N <- 400
  T <- 12
  
  set.seed(seed_seq[sim_number])
  # Generate the data Y
  Y <- matrix(rnorm(N*T, mean = alpha, sd = 1), nrow = N, ncol = T, byrow = FALSE)
  
  map <- map_partition(Y = Y, A_block = W, gamma_init = gamma_1, a1 = 1/T, a2 = T, nu_sigma = 3, lambda_sigma = 1, rho = 0.9)
  lam_1 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, gamma_init = gamma_1, a1 = 1/T, a2 = T, nu_sigma = 3, lambda_sigma = 1, rho = 0.9, lambda = 1)
  lam_10 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, gamma_init = gamma_1, a1 = 1/T, a2 = T, nu_sigma = 3, lambda_sigma = 1, rho = 0.9, lambda = 10)
  lam_100 <- ensm_cluster_mean(Y = Y, A_block = W, L = 10, gamma_init = gamma_1, a1 = 1/T, a2 = T, nu_sigma = 3, lambda_sigma = 1, rho = 0.9, lambda = 100)
  
  assign(paste0("map_sim", sim_number), map)
  assign(paste0("lam1_sim", sim_number), lam_1)
  assign(paste0("lam10_sim", sim_number), lam_10)
  assign(paste0("lam100_sim", sim_number), lam_100)
  
  save(list = paste0(c("map", "lam1", "lam10", "lam100"), "_sim", sim_number), file = paste0(dropbox_path, "/results/june12_sim", sim_number, ".RData"))
  
}

