//
//  local_search.cpp
//  Given an initial partition, form each island candidate and return the log-posterior
//  of each.
//  Created by Sameer Deshpande on 4/12/19.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include <vector>

using namespace std;

// [[Rcpp::export]]
Rcpp::List local_search(arma::vec ybar,
                        const int T,
                        const arma::mat A_block,
                        Rcpp::List gamma_init,
                        const double a1 = 1.0,
                        const double a2 = 1.0,
                        const double nu_sigma = 3,
                        const double lambda_sigma = 1,
                        const double rho = 0.99,
                        const double eta = 1.0,
                        const double eps = 1e-3)
{
  int n = ybar.size();
  Rcpp::Rcout << "n = " << n << endl;
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::Rcout << "Created gamma_0" << endl;
  
  split_info local_si;
  get_local(local_si, gamma_0, T, A_block, rho, a1, a2);
  
  int L = n+1;
  std::vector<LPPartition> particle_set(L);
  std::vector<double> w(L);
  particle_set[0] = new Partition(gamma_0);
  for(int i = 1; i < L; i++){
    particle_set[i] = new Partition(gamma_0);
    particle_set[i]->Split_Merge(local_si.split_k[i], local_si.new_clusters[i], local_si.nearest_neighbor[i], ybar, T, A_block, rho, a1, a2, eta);
    
    
    w[i] = 1/(n+1);
  }
  
  
  
  
  
  arma::vec tmp_log_post(n); // holds the log-posterior of the island candidate
  arma::vec tmp_log_like(n); // holds log-likelihood of each island candidate
  arma::vec tmp_log_post(n); // holds the log prior of each island candidate
  arma::uvec log_post_index(n);
  
  arma::vec log_post(n+1);
  arma::vec log_like(n+1);
  arma::vec log_prior(n+1);
  
  split_info local_si;
  
  double local_obj = 0.0;
  
  
  for(int i = 0; i < local_si.num_splits; i++){
    delete local_candidate;
    local_candidate = new Partition(gamma_0);
    local_candidate->Split_Merge(local_si.split_k[i], local_si.new_clusters[i], local_si.nearest_neighor[i], ybar, T, A_block, rho, a1, a2, eta);
    particle_set[i]->Copy_Partition(local_candidate);
  }
  
  
  
  
  log_post_index = arma::sort_index(tmp_log_post, "descend");
  
  log_post(0) = total_log_post(gamma_0, nu_sigma, lambda_sigma);
  log_like(0) = total_log_like(gamma_0, nu_sigma, lambda_sigma);
  log_prior(0) = total_log_prior(gamma_0);

  
  
  
}
