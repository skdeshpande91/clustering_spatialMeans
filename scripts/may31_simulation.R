library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/ensm_cluster_mean.cpp")
sourceCpp("src/map_cluster.cpp")


load("data/may31_sim_data.RData")




# Simulation 1

particle0_L5_lam10_sim1 <- ensm_cluster_mean(y1, A_block_large, L = 5, gamma_0_large, a1 = 1/10, a2 = 1/10, max_iter = 20, lambda = 10)
particle1_L5_lam10_sim1 <- ensm_cluster_mean(y1, A_block_large, L = 5, gamma_1_large, a1 = 1/10, a2 = 1/10, max_iter = 20, lambda = 10)

particle0_L5_lam100_sim1 <- ensm_cluster_mean(y1, A_block_large, L = 5, gamma_0_large, a1 = 1/10, a2 = 1/10, max_iter = 20, lambda = 100)
particle1_L5_lam100_sim1 <- ensm_cluster_mean(y1, A_block_large, L = 5, gamma_1_large, a1 = 1/10, a2 = 1/10, max_iter = 20, lambda = 100)



save(particle0_L5_lam10_sim1, particle1_L5_lam10_sim1,
     particle0_L5_lam100_sim1, particle1_L5_lam100_sim1,
     file = "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/results/may31_l5.RData")



particle1_L10_lam10_sim1 <- ensm_cluster_mean(y1, A_block_large, L = 10, gamma_1_large, 
                                            a1 = 1/10, a2 = 10, max_iter = 20, lambda = 10)




save(particle1_L5_lam1_sim1, file = "results/test.RData")