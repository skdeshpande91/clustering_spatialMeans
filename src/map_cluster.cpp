//
//  map_cluster.cpp
//  
//
//  Created by Sameer Deshpande on 4/7/19.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include <vector>

using namespace std;

// [[Rcpp::export]]
Rcpp::List map_partition(arma::vec ybar,
                         const int T,
                         const arma::mat A_block,
                         Rcpp::List gamma_init,
                         const double a1 = 1.0,
                         const double a2 = 1.0,
                         const double nu_sigma = 3,
                         const double lambda_sigma = 1,
                         const double rho = 0.99,
                         const double eta = 1.0,
                         const int max_iter = 10, const double eps = 1e-3,
                         const double split_frac = 0.1)
{
  int n = ybar.size();
  Rcpp::Rcout << "n = " << n << endl;
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::Rcout << "Created gamma_0" << endl;
  
  // for the main loop, we need the following quantities
  LPPartition spec_split_candidate = new Partition(gamma_0); // for spectral splits
  LPPartition tail_split_candidate = new Partition(gamma_0); // for tail splits
  LPPartition km_split_candidate = new Partition(gamma_0); // for km splits
  LPPartition merge_candidate = new Partition(gamma_0); // for merges
  LPPartition border_candidate = new Partition(gamma_0);
  LPPartition island_candidate = new Partition(gamma_0);
  
  split_info spec_si;
  split_info tail_si;
  split_info km_si;
  split_info bi;
  split_info isl_i;
  merge_info mi;
  
  double spec_split_obj = 0.0;
  double tail_split_obj = 0.0;
  double km_split_obj = 0.0;
  double merge_obj = 0.0;
  double border_obj = 0.0;
  double island_obj = 0.0;
  double accepted_obj = 0.0;
  
  int spec_split_flag = 1;
  int tail_split_flag = 1;
  int km_split_flag = 1;
  int merge_flag = 1;
  int border_flag = 1;
  int island_flag = 1;
  
  int iter = 0;
  double old_objective = 0.0;
  double objective = 0.0;
  int flag = 0;
  LPPartition gamma_hat = new Partition(gamma_0); // the partition we constantly update
  
  objective = total_log_post(gamma_hat, nu_sigma, lambda_sigma);
  while( (iter < max_iter) & (flag == 0)){
    Rcpp::Rcout << "Starting iter = " << iter << endl;
    old_objective = objective;
    Rcpp::checkUserInterrupt();
    spec_split_flag = 1;
    tail_split_flag = 1;
    km_split_flag = 1;
    border_flag = 1;
    merge_flag = 1;
    island_flag = 1;
    
    // Spectral Split
    get_spectral_split(spec_si, gamma_hat, T, A_block, rho, a1, a2, 1000);
    delete spec_split_candidate;
    spec_split_candidate = new Partition(gamma_hat);
    best_split_map(spec_si, spec_split_candidate, gamma_hat, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, true);
    spec_split_obj = total_log_post(spec_split_candidate, nu_sigma, lambda_sigma);
    spec_split_flag = Partition_Equal(spec_split_candidate, gamma_hat);
    accepted_obj = spec_split_obj;
    
    // Tail splits
    get_tail_split(tail_si, gamma_hat, T, A_block, rho, a1, a2, 0.025);
    delete tail_split_candidate;
    tail_split_candidate = new Partition(gamma_hat);
    best_split_map(tail_si, tail_split_candidate, gamma_hat, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, true);
    tail_split_obj = total_log_post(tail_split_candidate, nu_sigma, lambda_sigma);
    tail_split_flag = Partition_Equal(tail_split_candidate, gamma_hat);
    if(tail_split_obj > accepted_obj) accepted_obj = tail_split_obj;
    
    // KM splits
    get_km_split(km_si, gamma_hat, T, A_block, rho, a1, a2, 2500);
    delete km_split_candidate;
    km_split_candidate = new Partition(gamma_hat);
    best_split_map(km_si, km_split_candidate, gamma_hat, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, true);
    km_split_obj = total_log_post(km_split_candidate, nu_sigma, lambda_sigma);
    km_split_flag = Partition_Equal(km_split_candidate, gamma_hat);
    if(km_split_obj > accepted_obj) accepted_obj = km_split_obj;
    
    // Merges
    get_merge(mi, gamma_hat, A_block);
    delete merge_candidate;
    merge_candidate = new Partition(gamma_hat);
    best_merge_map(mi, merge_candidate, gamma_hat, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta);
    merge_obj = total_log_post(merge_candidate, nu_sigma, lambda_sigma);
    merge_flag = Partition_Equal(merge_candidate, gamma_hat);
    if(merge_obj > accepted_obj) accepted_obj = merge_obj;
    
    // Border
    get_border(bi, gamma_hat, T, A_block, rho, a1, a2);
    delete border_candidate;
    border_candidate = new Partition(gamma_hat);
    best_split_map(bi, border_candidate, gamma_hat, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, false);
    border_obj = total_log_post(border_candidate, nu_sigma, lambda_sigma);
    border_flag = Partition_Equal(border_candidate, gamma_hat);
    if(border_obj > accepted_obj) accepted_obj = border_obj;
    
    
    // Island
    get_island(isl_i, gamma_hat, T, A_block, rho, a2, a2, 0.05);
    delete island_candidate;
    island_candidate = new Partition(gamma_hat);
    best_split_map(isl_i, island_candidate, gamma_hat, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, false);
    island_obj = total_log_post(island_candidate, nu_sigma, lambda_sigma);
    island_flag = Partition_Equal(island_candidate, gamma_hat);
    if(island_obj > accepted_obj) accepted_obj = island_obj;
    
    if(accepted_obj == spec_split_obj){
      if(spec_split_flag == 0) gamma_hat->Copy_Partition(spec_split_candidate);
      else flag = 1; // best move is to not move at all
    } else if(accepted_obj == tail_split_obj){
      if(tail_split_flag == 0) gamma_hat->Copy_Partition(tail_split_candidate);
      else flag = 1;
    } else if(accepted_obj == km_split_obj){
      if(km_split_obj == 0) gamma_hat->Copy_Partition(km_split_candidate);
      else flag = 1;
    } else if(accepted_obj == merge_obj){
      if(merge_flag == 0) gamma_hat->Copy_Partition(merge_candidate);
      else flag = 1;
    } else if(accepted_obj == border_obj){
      if(border_flag == 0) gamma_hat->Copy_Partition(border_candidate);
      else flag = 1;
    } else if(accepted_obj == island_obj){
      if(island_flag == 0) gamma_hat->Copy_Partition(island_candidate);
      else flag = 1;
    }
    objective = total_log_post(gamma_hat, nu_sigma, lambda_sigma);
    Rcpp::Rcout << "objective = " << objective << "  old_objective = " << old_objective << " % diff = " << 100.0 * abs( (objective - old_objective)/objective) << endl;
    iter++;
  } // closes while loop
  
  std::vector<LPPartition> tmp_particle_set(1);
  tmp_particle_set[0] = gamma_hat;
  Rcpp::List gamma_out;
  format_particle_set(tmp_particle_set, gamma_out);
  Rcpp::List results;
  results["gamma_hat"] = gamma_out;
  
  return(results);
}
