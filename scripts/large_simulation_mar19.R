library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/ensm_cluster_mean2.cpp")

source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

load("data/scaled_alpha_large.RData")
load("data/partitions_large.RData")
load("data/large_example_1.RData")
load("data/large_example_2.RData")
load("data/large_example_3.RData")

lambda_seq <- c(400, 200, 100, 50, 25, 10, 5, 2, 1)
image_path <- "~/Dropbox/Particle EM Spatial Clustering/images/sameer_sim_mar20/"

# 5 particles, 
for(lambda in lambda_seq){
  tmp0 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large,
                            a1 = 1/10, a2 = 10, nu_sigma = 3, lambda_sigma = 1, lambda = lambda, L = 5)
  
  tmp1 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large,
                             a1 = 1/10, a2 = 10, nu_sigma = 3, lambda_sigma = 1, lambda = lambda, L = 5)
  
  png(paste0(image_path, "large2_L5_lam", lambda, "_0.png"), width = 6, height= 6, units= "in", res = 300)
  plot_particle_set_grid(tmp0, A_block_large, max_value = max_value)
  dev.off()
  
  png(paste0(image_path, "large2_L5_lam", lambda, "_1.png"), width = 6, height= 6, units= "in", res = 300)
  plot_particle_set_grid(tmp1, A_block_large, max_value = max_value)
  dev.off()
  
  assign(paste0("large2_L5_lam", lambda, "_0"), tmp0)
  assign(paste0("large2_L5_lam", lambda, "_1"), tmp1)
  
}

# 10 particles, 
for(lambda in lambda_seq){
  tmp0 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large,
                             a1 = 1/10, a2 = 10, nu_sigma = 3, lambda_sigma = 1, lambda = lambda, L = 10)
  
  tmp1 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large,
                             a1 = 1/10, a2 = 10, nu_sigma = 3, lambda_sigma = 1, lambda = lambda, L = 10)
  
  png(paste0(image_path, "large2_L10_lam", lambda, "_0.png"), width = 6, height= 6, units= "in", res = 300)
  plot_particle_set_grid(tmp0, A_block_large, max_value = max_value)
  dev.off()
  
  png(paste0(image_path, "large2_L10_lam", lambda, "_1.png"), width = 6, height= 6, units= "in", res = 300)
  plot_particle_set_grid(tmp1, A_block_large, max_value = max_value)
  dev.off()
  
  assign(paste0("large2_L10_lam", lambda, "_0"), tmp0)
  assign(paste0("large2_L10_lam", lambda, "_1"), tmp1)
  
}


# 25 particles
for(lambda in lambda_seq){
  tmp0 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large,
                             a1 = 1/10, a2 = 10, nu_sigma = 3, lambda_sigma = 1, lambda = lambda, L = 25)
  
  tmp1 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large,
                             a1 = 1/10, a2 = 10, nu_sigma = 3, lambda_sigma = 1, lambda = lambda, L = 25)
  
  png(paste0(image_path, "large2_L25_lam", lambda, "_0.png"), width = 6, height= 6, units= "in", res = 300)
  plot_particle_set_grid(tmp0, A_block_large, max_value = max_value)
  dev.off()
  
  png(paste0(image_path, "large2_L25_lam", lambda, "_1.png"), width = 6, height= 6, units= "in", res = 300)
  plot_particle_set_grid(tmp1, A_block_large, max_value = max_value)
  dev.off()
  
  assign(paste0("large2_L25_lam", lambda, "_0"), tmp0)
  assign(paste0("large2_L25_lam", lambda, "_1"), tmp1)
  
}

save.list <- c(paste0("large2_L5_lam", lambda_seq, "_0"), paste0("large2_L5_lam", lambda_seq, "_1"),
               paste0("large2_L10_lam", lambda_seq, "_0"), paste0("large2_L10_lam", lambda_seq, "_1"),
               paste0("large2_L25_lam", lambda_seq, "_0"), paste0("large2_L25_lam", lambda_seq, "_1"))

save(list = save.list, file = "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/march20_simulations.RData")


# Stuff from before


large2_L25_lam400_0 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large,
                                         a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 400, L = 25)
png("~/Dropbox/Particle EM Spatial Clustering/images/large2_L25_lam400_0.png", width = 6, height = 6, units = "in", res = 300)
plot_particle_set_grid(large2_L25_lam400_0, A_block_large, max_value = max_value)
dev.off()





large2_L5_lam100_0 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large,
                                         a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 100, L = 5)

png("~/Dropbox/Particle EM Spatial Clustering/images/large2_L25_lam400_0.png", width = 6, height = 6, units = "in", res = 300)
plot_particle_set_grid(large2_L25_lam400_0, A_block_large, max_value = max_value)
dev.off()

large2_L5_lam100_0 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_0_large,
                                         a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 100, L = 25)




large2_L5_lam400_1 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, gamma_init = gamma_1_large,
                                         a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, lambda = 400, L = 5)



large2_L25_lam10 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, L = 25, 
                                       gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                       lambda = 10)
large2_L25_lam100 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, L = 25, 
                                       gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                       lambda = 100)

png("~/Dropbox/Particle EM Spatial Clustering/images/large2_L25_lam10.png", width = 8, height = 8, units = "in", res = 300)
plot_particle_set_grid(large2_L25_lam10, A_block_large, max_value = max_value)
dev.off()


png("~/Dropbox/Particle EM Spatial Clustering/images/large2_L25_lam100.png", width = 8, height = 8, units = "in", res = 300)
plot_particle_set_grid(large2_L25_lam100, A_block_large, max_value = max_value)
dev.off()


test <-  ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, L = 5, 
                                       gamma_init = gamma_0_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                       lambda = 1)

large1_L25_lam1 <-  ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 25, 
                                       gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                       lambda = 1)
large2_L25_lam1 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, L = 25, 
                                      gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                      lambda = 1)
large3_L25_lam1 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, L= 25,
                                      gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                      lambda = 1)

large1_L25_lam10 <-  ensm_cluster_mean2(ybar_1_large[,1], T = 10, A_block = A_block_large, L = 25, 
                                       gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                       lambda = 10)

large3_L25_lam10 <- ensm_cluster_mean2(ybar_3_large[,1], T = 10, A_block = A_block_large, L= 25,
                                      gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                      lambda = 10)


large2_L25_lam10 <- ensm_cluster_mean2(ybar_2_large[,1], T = 10, A_block = A_block_large, L = 25, 
                                       gamma_init = gamma_1_large, a1 = 1/10, a2 = 10, nu_sigma = 1, lambda_sigma = 1, 
                                       lambda = 100)
