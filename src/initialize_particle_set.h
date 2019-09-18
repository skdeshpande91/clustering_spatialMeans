//
//  initialize_particle_set.hpp
//  
//
//  Created by Sameer Deshpande on 9/9/19.
//

#ifndef GUARD_initialize_particle_set_h
#define GUARD_initialize_particle_set_h

#include<RcppArmadillo.h>
#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"
#include "rng.h"
#include <stdio.h>

// init_id == 0: make everything equal to partition with everyone in the same cluster
// init_id == 1: make L/log(N) copies of each output of K-means
// init_id == 2: importance sample

void initialize_particle_set(std::vector<LPPartition> &particle_set, const int L, const arma::vec ybar, const double total_ss, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double eta, const double nu_sigma, const double lambda_sigma, const int init_id, RNG &gen);


#endif /* initialize_particle_set_hpp */
