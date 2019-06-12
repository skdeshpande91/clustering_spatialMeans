# Partition Functions
library(scales)
library(RColorBrewer)
col_list <- rev(brewer.pal(n = 5, name = "RdBu"))

sorted_equal <- function(v1, v2){
  return(identical(sort(as.numeric(v1)), sort(as.numeric(v2))))
}

# A function to check whether x is inside vector v. 
# we will use this with apply to check whether

index_in <- function(x, v){
  return(x %in% v)
}

# Configuration of the parttion
partition_config <- function(gamma, A_block){
  tmp <- sapply(gamma, FUN = length)
  singletons <- which(tmp == 1)
  non_singletons <- which(tmp != 1)
  K <- length(gamma) # number of clusters
  n <- nrow(A_block) # number of block groups
  n_singletons <- length(singletons)
  n_non_singletons <- length(non_singletons)
  Z <- matrix(0, nrow = n, ncol = n) # pairwise co-allocation matrix
  
  for(k in 1:K){
    Z[gamma[[k]], gamma[[k]]] <- 1
  }
  
  if(K > 1){ # identify which clusters are spatially adjacent
    A_cluster <- matrix(0, nrow = K, ncol = K)
    for(k in 1:(K-1)){
      for(kk in (k+1):K){
        A_cluster[k,kk] <- 1*any(A_block[gamma[[k]], gamma[[kk]]] == 1)
        A_cluster[kk,k] <- A_cluster[k,kk]
      }
    }
  } else{
    A_cluster <- NULL
  }
  return(list(K = K, config = tmp, 
              singletons = singletons, non_singletons = non_singletons,
              A_cluster = A_cluster, Z = Z))
}

# Ensures that the clusters are spatially connected
partition_modify <- function(gamma, A_block){
  gamma_new <- gamma
  K <- length(gamma)
  count <- K + 1
  
  to_be_removed <- c()
  for(k in 1:K){
    if(length(gamma[[k]])>1){
      cl <- gamma[[k]]
      g <- graph_from_adjacency_matrix(A_block[cl,cl]) # adjacency matrix of the cluster
      tmp <- components(g)$membership # finds the connected components of the cluster
      if(length(unique(tmp))>1){ # there are more than connected components
        for(x in unique(tmp)){ # loops over the component names
          gamma_new[[count]] <- cl[which(tmp == x)] # forms a new cluster out the connected components. adds to end of gamma
          count <- count + 1
        }
        to_be_removed <- c(to_be_removed, k) # will eventually delete the original cluster j
      }
    }
  }
  if(length(to_be_removed) > 0) gamma_new <- gamma_new[-to_be_removed]
  return(gamma_new)
}

partition_equal <- function(gamma1, gamma2, A_block){
  tmp_1 <- partition_config(gamma1, A_block)
  tmp_2 <- partition_config(gamma2, A_block)
  
  K_1 <- tmp_1$K
  K_2 <- tmp_2$K
  
  config_1 <- sort(tmp_1$config)
  config_2 <- sort(tmp_2$config)
  
  if(! identical(config_1, config_2)){ # if the cluster configurations are different, stop
    flag <- FALSE
  } else { # now check individual partition elements
    flag <- TRUE
    i <- 0
    while(i < K_1 & flag == TRUE){
      i <- i+1
      flag <- any(sapply(gamma2, FUN = sorted_equal, gamma1[[i]]))
    }
  }
  return(flag)
}

vi_distance <- function(N, gamma1, gamma2){
  # Create the N matrix
  k1 <- length(gamma1)
  k2 <- length(gamma2)
  counts <- matrix(0, nrow = k1, ncol = k2)
  for(k in 1:k1){
    for(kk in 1:k2){
      counts[k,kk] <- length(intersect(gamma1[[k]], gamma2[[kk]]))
    }
  }
  row_sums <- rowSums(counts)
  col_sums <- colSums(counts)
  dist <- sum(row_sums/N * log(row_sums/N)) + sum(col_sums/N * log(col_sums/N))
  for(k in 1:k1){
    for(kk in 1:k2){
      if(counts[k,kk] > 0){
        dist <- dist -2 * counts[k,kk]/N * log(counts[k,kk]/N)
      }
    }
  }
  
  return(dist)
}


binder_loss <- function(N, gamma1, gamma2){
  
  k1 <- length(gamma1)
  k2 <- length(gamma2)
  counts <- matrix(nrow = k1, ncol = k2)
  # need matrix of counts n_ij counting the number of indices that belong to cluster k in gamma1 and cluster kk in gamma2
  for(k in 1:k1){
    for(kk in 1:k2){
      counts[k,kk] <- length(intersect(gamma1[[k]], gamma2[[kk]]))
    }
  }
  row_sums <- rowSums(counts)
  col_sums <- colSums(counts)
  dist <- sum( (row_sums/N)^2) + sum( (col_sums/N)^2)
  for(k in 1:k1){
    for(kk in 1:k2){
      if(counts[k,kk] > 0){
        dist <- dist - 2 * (counts[k,kk]/N)^2
      }
    }
  }
  
  return(dist)
  
}
