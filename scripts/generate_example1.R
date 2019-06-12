# Generate data for Example 1
# We will generate data from a misspecified model

library(igraph)
library(MASS)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/ensm_cluster_mean.cpp")
sourceCpp("src/map_cluster.cpp")


source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

##############
# A small grid graph example
##############


n_side_large <- 20
g_large <- make_lattice(length = n_side_large, dim = 2)

A_block_large <- as_adj(g_large, type = "both")
A_block_large <- matrix(A_block_large, nrow = n_side_large^2)

N_large <- dim(A_block_large)[1]
w_sym_large <- A_block_large

############
# Divide grid into quadrants
############


quad_1_large <- rep(11:20, times = 10) + 20*rep(10:19, each = 10)
quad_2_large <- rep(1:10, times = 10) + 20*rep(10:19, each = 10)
quad_3_large <- rep(1:10, times = 10) + 20*rep(0:9, each = 10)
quad_4_large <- rep(11:20, times = 10) + 20*rep(0:9, each = 10)

# For the large grid, we need to define the actual clusters
cluster_1_large <- c(165, 166, 184:187, 204:207, 225, 226)
cluster_2_large <- c(quad_2_large, quad_3_large)
cluster_2_large <- cluster_2_large[!cluster_2_large %in% cluster_1_large]
cluster_3_large <- quad_4_large
cluster_4_large <- quad_1_large


# Set up 5 partitions for the large grid case
gamma_0_large <- list(cluster_1_large, cluster_2_large, cluster_3_large, cluster_4_large)
gamma_1_large <- list(c(quad_1_large, quad_2_large, quad_3_large, quad_4_large))
gamma_2_large <- list(c(quad_1_large, quad_2_large), c(quad_3_large, quad_4_large))
gamma_3_large <- list(quad_1_large, quad_2_large, c(quad_3_large, quad_4_large))
gamma_4_large <- list(quad_1_large, quad_2_large, quad_3_large, quad_4_large)

# Generate the alpha's and y
tmp_large <- partition_config(gamma_0_large, A_block_large)
K_large <- tmp_large$K
config_large <- tmp_large$config

set.seed(129)
sigma2 <- 1/rgamma(1, shape = 3/2, rate = 3/2) # sigma2 ~ 3/chi-square_3

alpha_bar1 <- c(0, -5, 3, 1)
alpha_bar2 <- c(0, -2, 3, 1)

#alpha_bar <- c(0, 3, -2, 1.25)
#alpha_bar <- c(0, -2, 3, 1)
T <- 12
alpha1 <- rep(NA, times = N_large)
alpha2 <- rep(NA, times = N_large)

y1 <- matrix(NA, nrow = N_large, ncol = T)
y2 <- matrix(NA, nrow = N_large, ncol = T)

set.seed(1991)
for(k in 1:K_large){
  alpha1[gamma_0_large[[k]]] <- runif(config_large[k], alpha_bar1[k] - 0.75, alpha_bar1[k] + 0.75)
  alpha2[gamma_0_large[[k]]] <- runif(config_large[k], alpha_bar2[k] - 0.75, alpha_bar2[k] + 0.75)
  for(i in gamma_0_large[[k]]){
    y1[i,] <- rnorm(T, mean = alpha1[i], sd = sqrt(sigma2))
    y2[i,] <- rnorm(T, mean = alpha2[i], sd = sqrt(sigma2))
  }
}

#max_value <- max(c(abs(y1),abs(alpha1), abs(y2), abs(alpha2)))
max_value <- max(abs(alpha1), abs(alpha2), abs(rowMeans(y1)), abs(rowMeans(y2)))


save(gamma_0_large, gamma_1_large, gamma_2_large, gamma_3_large, gamma_4_large,
     A_block_large, alpha_bar1, alpha_bar2, alpha1, alpha2, y1, y2, sigma2, max_value,
     file = "data/may31_sim_data.RData")



plot_partition_grid(gamma_1_large, A_block_large, values = alpha1, max_value = max_value)
plot_partition_grid(gamma_1_large, A_block_large, values = alpha2, max_value = max_value)

test1_map0 <- map_partition(y1, A_block_large, gamma_0_large, a1 = 1/10, a2 = 10, max_iter = 20)
test_map1 <- map_partition(y, A_block_large, gamma_1_large, a1 = 1/10, a2 = 10, max_iter = 20)

test_particle0 <- ensm_cluster_mean(y, A_block_large, L = 5, gamma_0_large,
                                    a1 = 1/10, a2 = 10, max_iter = 20)

test_particle1 <- ensm_cluster_mean(y, A_block_large, L = 5, gamma_1_large,
                                    a1 = 1/10, a2 = 10, max_iter = 20)

test_particle0_lam100 <- ensm_cluster_mean(y, A_block_large, L = 5, gamma_0_large,
                                           a1 = 1/10, a2 = 10, max_iter= 20, 
                                           lambda = 100)

test_particle1_lam100 <- ensm_cluster_mean(y, A_block_large, L = 5, gamma_1_large,
                                            a1 = 1/10, a2 = 10, max_iter= 20, 
                                           lambda = 100)
test1_L10_lam100 <- ensm_cluster_mean(y, A_block_large, L = 10, gamma_1_large, 
                                      a1 = 1/10, a2 = 10, max_iter = 25, lambda = 100)
# objective starting from gamma_0_large; -16481.2
# 
path <- "~/Dropbox/Particle EM Spatial Clustering/sameer_sim/spring2019/"

save(gamma_0_large, gamma_1_large, gamma_2_large, gamma_3_large, gamma_4_large,
     A_block_large, alpha_bar, alpha, y, sigma2, max_value,
     test_map0, test_map1, test_particle0_lam1, test_particle0_lam100,
     test_particle1_lam1, test_particle1_lam100, test1_L10_lam100,
     file = paste0(path, "/may16_example1.RData"))

