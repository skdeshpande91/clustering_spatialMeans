//
//  test_init_particles.cpp
//  
//
//  Created by Sameer Deshpande on 9/9/19.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include "initialize_particle_set.h"
#include <vector>
#include <ctime>


// [[Rcpp::export]]
Rcpp::List test_init_particles(arma::mat Y,
                               const arma::mat A_block,
                               const int L = 10,
                               const double a1 = 1.0,
                               const double a2 = 1.0,
                               const double nu_sigma = 3,
                               const double lambda_sigma = 1,
                               const double rho = 0.99,
                               const double lambda = 1.0, const double eta = 1.0,
                               const int max_iter = 10, const double eps = 1e-3,
                               bool verbose = false)
{
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }
  
  std::vector<LPPartition> particle_set(L);
  std::vector<double> w(L);
  for(int l = 0; l < L; l++) w[l] = 1.0/ ( (double) L);
  
  Rcpp::Rcout << "About to start initialize_particle_set" << std::endl;
  for(int l = 0; l < L; l++) particle_set[l] = new Partition();
  initialize_particle_set(particle_set, L, ybar, total_ss, T, A_block, rho,  a1, a2, eta, nu_sigma, lambda_sigma);
  
  // Find the unique particles
  std::vector<std::vector<int> > particle_map;
  std::vector<double> pstar;
  std::vector<int> counts;
  get_unik_particles(particle_map, pstar, counts, particle_set, w);

  int L_star = particle_map.size();
  std::vector<LPPartition> unik_particles(L_star);
  std::vector<double> log_like(L_star);
  std::vector<double> log_prior(L_star);
  std::vector<double> log_post(L_star);
  for(int l = 0; l < L_star; l++){
    unik_particles[l] = new Partition(particle_set[particle_map[l][0]]);
    log_like[l] = total_log_like(unik_particles[l], total_ss, T, nu_sigma, lambda_sigma);
    log_prior[l] = total_log_prior(unik_particles[l]);
    log_post[l] = total_log_post(unik_particles[l], total_ss, T, nu_sigma, lambda_sigma);
  }
  
  
  Rcpp::List unik_particles_out;
  format_particle_set(unik_particles, unik_particles_out);

  
  Rcpp::List results;
  //results["ybar"] = ybar;
  results["particles"] = unik_particles_out;
  results["pstar"] = pstar;
  results["counts"] = counts;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  return(results);

  
  
}
