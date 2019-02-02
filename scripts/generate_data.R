#################
# Generate two sets of data to test the code
#################


library(igraph)
library(MASS)
source("scripts/partition_functions.R")
source("scripts/plot_partition.R")

##############
# A small grid graph example
##############

n_side_small <- 10
n_side_large <- 20
g_small <-make_lattice(length = n_side_small, dim = 2)
g_large <- make_lattice(length = n_side_large, dim = 2)

A_block_small <- as_adj(g_small, type = "both")
A_block_large <- as_adj(g_large, type = "both")
A_block_small <- matrix(A_block_small, nrow = n_side_small^2)
A_block_large <- matrix(A_block_large, nrow = n_side_large^2)

N_small <- dim(A_block_small)[1] # number of block-groups
N_large <- dim(A_block_large)[1]
w_sym_small <- A_block_small
w_sym_large <- A_block_large

############
# Divide grid into quadrants
############

quad_1_small <- rep(6:10, times = 5) + 10*rep(5:9, each = 5)
quad_2_small <- rep(1:5, times = 5) + 10*rep(5:9, each = 5)
quad_3_small <- rep(1:5, times = 5) + 10*rep(0:4, each = 5)
quad_4_small <- rep(6:10, times = 5) + 10*rep(0:4, each = 5)

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

# Set up 5 partitions for the small grid case
gamma_0_small <- list(c(quad_2_small, quad_3_small), quad_4_small, quad_1_small)
gamma_1_small <- list(c(quad_1_small, quad_2_small, quad_3_small, quad_4_small))
gamma_2_small <- list(c(quad_1_small, quad_2_small), c(quad_3_small, quad_4_small))
gamma_3_small <- list(quad_1_small, quad_2_small, c(quad_3_small, quad_4_small))
gamma_4_small <- list(quad_1_small, quad_2_small, quad_3_small, quad_4_small)

# Set up 5 partitions for the large grid case
gamma_0_large <- list(cluster_1_large, cluster_2_large, cluster_3_large, cluster_4_large)
gamma_1_large <- list(c(quad_1_large, quad_2_large, quad_3_large, quad_4_large))
gamma_2_large <- list(c(quad_1_large, quad_2_large), c(quad_3_large, quad_4_large))
gamma_3_large <- list(quad_1_large, quad_2_large, c(quad_3_large, quad_4_large))
gamma_4_large <- list(quad_1_large, quad_2_large, quad_3_large, quad_4_large)

save(gamma_0_small, gamma_1_small, gamma_2_small, gamma_3_small, gamma_4_small, A_block_small, file = "data/partitions_small.RData")
save(gamma_0_large, gamma_1_large, gamma_2_large, gamma_3_large, gamma_4_large, A_block_large, file = "data/partitions_large.RData")
########################
# Generate alpha
########################

# For small grid examples
rho <- 0.99
a1_small <- 1/10 # variance within cluster
a2_small <-  10 # cluster mean variance
a1_large <- 1/20
a2_large <- 20

tmp_small <- partition_config(gamma_0_small, A_block_small)
tmp_large <- partition_config(gamma_0_large, A_block_large)
K_small <- tmp_small$K
K_large <- tmp_large$K
config_small <- tmp_small$config
config_large <- tmp_large$config

alpha_bar_1_small <- c(-5, 5, 0)
alpha_bar_2_small <- c(-2, 2, 0)
alpha_bar_3_small <- c(-5, 2, -1)

alpha_bar_1_large <- c(0, -5, 5,0)
alpha_bar_2_large <- c(0, -2, 2, 0)
alpha_bar_3_large <- c(-2, -5, 5, 1)


# We will work with data from 2006 to 2015
T <- 10
n_sim <- 500
seed_list <- c(12991, 72460, 60357)
for(sim_number in 1:3){
  alpha_bar_small <- get(paste0("alpha_bar_", sim_number, "_small"))
  alpha_bar_large <- get(paste0("alpha_bar_", sim_number, "_large"))
  set.seed(seed_list[sim_number])
  #sigma2 <- 1/rgamma(K, 3/2, 1/2) # sigma2 ~ inv. chi-square_3
  #sigma2 <- 1/rgamma(1, 1, 1/2)
  sigma2 <- 1/rgamma(1, 3/2, 3/2) # sigma2 ~ 3/chi-square_3
  alpha_small <- rep(NA, times = N_small)
  alpha_large <- rep(NA, times = N_large)
  
  ybar_small <- matrix(NA, nrow = N_small, ncol = n_sim)
  ybar_large <- matrix(NA, nrow = N_large, ncol = n_sim)
  

  
  for(k in 1:K_small){
    n_k <- config_small[k]
    A_block_k <- A_block_small[gamma_0_small[[k]], gamma_0_small[[k]]] # adjacency matrix of the cluster
    D <- diag(rowSums(A_block_k))
    A_star_k <- D - A_block_k # unweighted laplacian 
    Omega_alpha <- rho * A_star_k + (1 - rho) * diag(n_k)
    Sigma_alpha <- solve(Omega_alpha)
    set.seed(seed_list[sim_number] + k)
    alpha_small[gamma_0_small[[k]]] <- mvrnorm(n = 1, mu = rep(alpha_bar_small[k], times = n_k), Sigma = a1_small * sigma2 * Sigma_alpha)
    ybar_small[gamma_0_small[[k]],] <- rnorm(n = n_k * n_sim, mean = alpha_small[gamma_0_small[[k]]], sd = sqrt(sigma2)/T)
  }
  for(k in 1:K_large){
    n_k <- config_large[k]
    A_block_k <- A_block_large[gamma_0_large[[k]], gamma_0_large[[k]]]
    D <- diag(rowSums(A_block_k))
    A_star_k <- D - A_block_k
    Omega_alpha <- rho*A_star_k + (1-rho)*diag(n_k)
    Sigma_alpha <- solve(Omega_alpha)
    set.seed(seed_list[sim_number] + k)
    alpha_large[gamma_0_large[[k]]] <- mvrnorm(n = 1, mu = rep(alpha_bar_large[k], times = n_k), Sigma = a1_large * sigma2 * Sigma_alpha)
    ybar_large[gamma_0_large[[k]],] <- rnorm(n = n_k * n_sim, mean = alpha_large[gamma_0_large[[k]]], sd = sqrt(sigma2)/T)
  }
  
  assign(paste0("alpha_", sim_number, "_small"), alpha_small)
  assign(paste0("alpha_", sim_number, "_large"), alpha_large)
  assign(paste0("ybar_", sim_number, "_small"), ybar_small)
  assign(paste0("ybar_", sim_number, "_large"), ybar_large)
  
  assign(paste0("sigma2_", sim_number), sigma2)
}

# Save all of the alphas and ybar's in .RData files. These will NOT sit on GitHub
save(alpha_bar_1_small, alpha_1_small, ybar_1_small, sigma2_1, file = "data/small_example_1.RData")
save(alpha_bar_2_small, alpha_2_small, ybar_2_small, sigma2_2, file = "data/small_example_2.RData")
save(alpha_bar_3_small, alpha_3_small, ybar_3_small, sigma2_3, file = "data/small_example_3.RData")

save(alpha_bar_1_large, alpha_1_large, ybar_1_large, sigma2_1, file = "data/large_example_1.RData")
save(alpha_bar_2_large, alpha_2_large, ybar_2_large, sigma2_2, file = "data/large_example_2.RData")
save(alpha_bar_3_large, alpha_3_large, ybar_3_large, sigma2_3, file = "data/large_example_3.RData")

# Save some alpha values as a csv ... these will sit on GitHub
alpha_small_df <- data.frame("sim1" = alpha_1_small, "sim2" = alpha_2_small, "sim3" = alpha_3_small)
alpha_large_df <- data.frame("sim1" = alpha_1_large, "sim2" = alpha_2_large, "sim3" = alpha_3_large)

y_small_df <- data.frame("sim1" = ybar_1_small[,1], "sim2" = ybar_2_small[,1], "sim3" = ybar_3_small[,1])
y_large_df <- data.frame("sim1" = ybar_1_large[,1], "sim2" = ybar_2_large[,1], "sim3" = ybar_3_large[,1])

write.csv(alpha_small_df, col.names = TRUE, row.names = FALSE, file = "data/alpha_small.csv")
write.csv(alpha_large_df, col.names = TRUE, row.names = FALSE, file = "data/alpha_large.csv")
write.csv(y_small_df, col.names = TRUE, row.names = FALSE, file = "data/ybar_small.csv")
write.csv(y_large_df, col.names = TRUE, row.names = FALSE, file = "data/ybar_large.csv")

write.csv(A_block_small, col.names = FALSE, row.names = FALSE, file = "data/A_block_small.csv")
write.csv(A_block_large, col.names = FALSE, row.names = FALSE, file = "data/A_block_large.csv")

# Get the correct scaling: this is primarily for plotting purposes
max_value_small <- ceiling(max(c(abs(ybar_1_small), abs(ybar_2_small), abs(ybar_3_small), 
                           abs(alpha_1_small), abs(alpha_2_small), abs(alpha_3_small))))

max_value_large <- ceiling(max(c(abs(ybar_1_large), abs(ybar_2_large), abs(ybar_3_large), 
                                 abs(alpha_1_large), abs(alpha_2_large), abs(alpha_3_large) )))

max_value <- max(max_value_small, max_value_large)
scaled_alpha_1_small <- rescale(alpha_1_small, to = c(0,1), from = c(-1*max_value, 1*max_value))
scaled_alpha_2_small <- rescale(alpha_2_small, to = c(0,1), from = c(-1*max_value, 1*max_value))
scaled_alpha_3_small <- rescale(alpha_3_small, to = c(0,1), from = c(-1*max_value, 1*max_value))

scaled_alpha_1_large <- rescale(alpha_1_large, to = c(0,1), from = c(-1*max_value, 1*max_value))
scaled_alpha_2_large <- rescale(alpha_2_large, to = c(0,1), from = c(-1*max_value, 1*max_value))
scaled_alpha_3_large <- rescale(alpha_3_large, to = c(0,1), from = c(-1*max_value, 1*max_value))

save(scaled_alpha_1_small, scaled_alpha_2_small, scaled_alpha_3_small, max_value, file = "data/scaled_alpha_small.RData")
save(scaled_alpha_1_large, scaled_alpha_2_large, scaled_alpha_3_large, max_value, file = "data/scaled_alpha_large.RData")
