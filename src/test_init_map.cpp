//
//  test_init_map.cpp
//  
//
//  Created by Sameer Deshpande on 9/23/19.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include "initialize_particle_set.h"
#include "rng.h"
#include <vector>
#include <ctime>

using namespace std;

// [[Rcpp::export]]
Rcpp::List test_init_map(arma::mat Y,
                         const arma::mat A_block,
                         const double a1 = 1.0,
                         const double a2 = 1.0,
                         const double nu_sigma = 3,
                         const double lambda_sigma = 1,
                         const double rho = 0.99,
                         const double eta = 1.0,
                         const int max_iter = 10, const double eps = 1e-3,
                         const double split_frac = 0.1,
                         bool verbose = false)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }
  LPPartition gamma_init = new Partition();
  initialize_particle(gamma_init, ybar, total_ss, T, A_block, rho, a1, a2, eta, nu_sigma, lambda_sigma, gen);
  gamma_init->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
  
  
  Rcpp::List results;
  results["Y"] = Y;
  return(results);
}
