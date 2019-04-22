//
//  ensm_cluster_mean.cpp
//
//  Updated by Sameer Deshpande on 1 April 2019.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include <vector>
#include <ctime>

using namespace arma;
using namespace std;


using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List ensm_cluster_mean(arma::mat Y,
                              const arma::mat A_block,
                              const int L,
                              Rcpp::List gamma_init,
                              const double a1 = 1.0,
                              const double a2 = 1.0,
                              const double nu_sigma = 3,
                              const double lambda_sigma = 1,
                              const double rho = 0.99,
                              const double lambda = 1.0, const double eta = 1.0,
                              const int max_iter = 10, const double eps = 1e-3,
                              const double split_frac = 0.1)
{
  
  int n = Y.n_rows;
  int T = Y.n_cols;
  
  arma::vec ybar(n);
  double total_ss = 0;
  for(int i = 0; i < n; i++){
    ybar(i) = arma::mean(Y.row(i));
    total_ss += (T-1) * arma::var(Y.row(i));
  }

  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::Rcout << "Created gamma_0" << endl;
  
  // create a particle set
  std::vector<LPPartition> particle_set(L);
  std::vector<double> w(L);
  for(int l = 0; l < L; l++){
    particle_set[l] = new Partition(gamma_0);
    w[l] = 1.0/( (double) L);
  }
  
  
  Rcpp::Rcout << "Initialized the particle set" << endl;
  
  // now that we are tracking the trajectory, we no longer need this
  /*
  std::vector<std::vector<int> > init_particle_map;
  std::vector<double> init_pstar;
  std::vector<int> init_counts;
  get_unik_particles(init_particle_map, init_pstar, init_counts, particle_set, w);
  
  int init_L_star = init_particle_map.size();
  std::vector<double> init_log_like(init_L_star);
  std::vector<double> init_log_prior(init_L_star);
  std::vector<double> init_log_post(init_L_star);
  
  std::vector<LPPartition> init_unik_particles(init_L_star);
  for(int l = 0; l < init_L_star; l++){
    init_unik_particles[l] = new Partition(particle_set[init_particle_map[l][0]]);
    init_log_like[l] = total_log_like(init_unik_particles[l], nu_sigma, lambda_sigma);
    init_log_prior[l] = total_log_prior(init_unik_particles[l]);
    init_log_post[l] = total_log_post(init_unik_particles[l], nu_sigma, lambda_sigma);
  }
  */
  
  // for the main loop we need the following quantitites
  LPPartition spec_split_candidate = new Partition(gamma_0); // for the spectral splits
  LPPartition tail_split_candidate = new Partition(gamma_0); // for the tail splits
  LPPartition km_split_candidate = new Partition(gamma_0); // for km splits
  LPPartition merge_candidate = new Partition(gamma_0); // for the merge candidates
  LPPartition border_candidate = new Partition(gamma_0); // for the border moves
  LPPartition island_candidate = new Partition(gamma_0); // for the island candidates
  LPPartition local_candidate = new Partition(gamma_0); // for local candidate
  
  split_info spec_si; // holds the information for spectral splits
  split_info tail_si; // holds the information for the tail splits
  split_info km_si; // holds the information for the km splits
  split_info bi ; // holds information for border moves
  split_info isl_i; // holds information for island moves
  merge_info mi; // holds informatino for the merge moves
  split_info local_si; // holds information for the local moves
  
  double spec_split_obj = 0.0;
  double tail_split_obj = 0.0;
  double km_split_obj = 0.0;
  double merge_obj = 0.0;
  double border_obj = 0.0;
  double island_obj = 0.0;
  double local_obj = 0.0;
  double accepted_obj = 0.0;
  
  int spec_split_flag = 1; // essentially will hold the value of Partition_Equal(spec_split_candidate, particle_set[l])
  int tail_split_flag = 1; // holds value of Partition_Equal(tail_split_candidate, particle_set[l])
  int km_split_flag = 1; // holds value of Partition_Equal(km_split_candidate, particle_set[l])
  int merge_flag = 1; // holds value of Partition_Equal(merge_candidate, particle_set[l]);
  int border_flag = 1; // holds value of Partition_Equal(border_candidate, particle_set[l]);
  int island_flag = 1; // holds value of Partition_Equal(island_candidate, particle_set[l])
  int local_flag = 1; // holds value of Partition_Equal(local_candidate, particle_set[l]);
  
  bool try_local = true;
  
  int iter = 0;
  int conv_counter = 0; // count the number of particles that remain the same
  int flag = 0; //
  
  double old_objective = 0.0;
  double objective = 0.0;
  objective = lambda * Entropy(0, particle_set[0], particle_set, w);
  //for(int l = 0; l < L; l++) objective += w[l] * total_log_post(particle_set[l], nu_sigma, lambda_sigma);
  for(int l = 0; l < L; l++) objective += w[l] * total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
  
  // Prepare containers to track the trajectories of the particle system
  
  std::vector<double> objective_trajectory;
  objective_trajectory.push_back(objective);
  
  std::vector<Rcpp::List> particle_set_trajectory;
  Rcpp::List tmp_list;
  format_particle_set(particle_set, tmp_list);
  particle_set_trajectory.push_back(tmp_list);
  
  std::vector<double> tmp_log_like(L); // holds
  std::vector<double> tmp_log_prior(L);
  std::vector<double> tmp_log_post(L);
  std::vector<std::vector<double> > tmp_alpha_hat(L, std::vector<double>(n)); // holds estimate of alpha in each iteration
  
  std::vector<std::vector<double> > log_like_trajectory;
  std::vector<std::vector<double> > log_prior_trajectory;
  std::vector<std::vector<double> > log_post_trajectory;
  std::vector<std::vector<std::vector<double> > > alpha_trajectory;
  
  for(int l = 0; l < L; l++){
    tmp_log_like[l] = total_log_like(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
    tmp_log_prior[l] = total_log_prior(particle_set[l]);
    tmp_log_post[l] = total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
    for(int i = 0; i < n; i++) tmp_alpha_hat[l][i] = particle_set[l]->alpha_hat[i];
  }
  
  log_like_trajectory.push_back(tmp_log_like);
  log_prior_trajectory.push_back(tmp_log_prior);
  log_post_trajectory.push_back(tmp_log_post);
  alpha_trajectory.push_back(tmp_alpha_hat);

  time_t tp;
  int time1 = time(&tp);
  
  while((iter < max_iter) & (flag == 0)){
    Rcpp::Rcout << "[ensm_cluster_mean]: Starting iter " << iter << endl;
    old_objective = objective;
    
    conv_counter = 0; // counts the number of particles unchanged in our sweep
    for(int l = 0; l < L; l++){
      Rcpp::Rcout << "Starting to update particle " << l << endl;
      Rcpp::checkUserInterrupt();
      spec_split_flag = 1;
      tail_split_flag = 1;
      km_split_flag = 1;
      border_flag = 1;
      merge_flag = 1;
      island_flag = 1;
      local_flag = 1;
      try_local = true;
      
      // Spectral Splits
      get_spectral_split(spec_si, particle_set[l], T, A_block, rho, a1, a2, 500);
      delete spec_split_candidate;
      spec_split_candidate = new Partition(particle_set[l]);
      best_split(spec_si, spec_split_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda, true);
      
      spec_split_obj = w[l]*total_log_post(spec_split_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda * Entropy(l, spec_split_candidate, particle_set, w);
      spec_split_flag = Partition_Equal(spec_split_candidate, particle_set[l]);
      if(spec_split_flag == 0) try_local = false; // we might accept a spectral split so don't try all of the local moves
      accepted_obj = spec_split_obj;
      
      // Tail splits
      get_tail_split(tail_si, particle_set[l], T, A_block, rho, a1, a2, 0.025);
      delete tail_split_candidate;
      tail_split_candidate = new Partition(particle_set[l]);
      best_split(tail_si, tail_split_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda, true);
      
      tail_split_obj = w[l]*total_log_post(tail_split_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda * Entropy(l, tail_split_candidate, particle_set,w);
      tail_split_flag = Partition_Equal(tail_split_candidate, particle_set[l]);
      if(tail_split_flag == 0) try_local = false;
      if(tail_split_obj > accepted_obj) accepted_obj = tail_split_obj;
      
      // KM splits
      get_km_split(km_si, particle_set[l], T, A_block, rho, a1, a2, 1000);
      delete km_split_candidate;
      km_split_candidate = new Partition(particle_set[l]);
      best_split(km_si, km_split_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda, true);
      
      km_split_obj = w[l]*total_log_post(km_split_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(l, km_split_candidate, particle_set,w);
      km_split_flag = Partition_Equal(km_split_candidate, particle_set[l]);
      if(km_split_flag == 0) try_local = false;
      if(km_split_obj > accepted_obj) accepted_obj = km_split_obj;
      // Merges
      get_merge(mi, particle_set[l], A_block);
      delete merge_candidate;
      merge_candidate = new Partition(particle_set[l]);
      best_merge(mi, merge_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      merge_obj = w[l]*total_log_post(merge_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(l, merge_candidate, particle_set, w);
      merge_flag = Partition_Equal(merge_candidate, particle_set[l]);
      if(merge_flag == 0) try_local = false;
      if(merge_obj > accepted_obj) accepted_obj = merge_obj;
      
      // Border
      get_border(bi, particle_set[l], T, A_block, rho, a1, a2);
      delete border_candidate;
      border_candidate = new Partition(particle_set[l]);
      best_split(bi, border_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda, false);
      border_obj = w[l]*total_log_post(border_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(l, border_candidate, particle_set,w);
      border_flag = Partition_Equal(border_candidate, particle_set[l]);
      if(border_flag == 0) try_local = false;
      if(border_obj > accepted_obj) accepted_obj = border_obj;
      
      // Island
      get_island(isl_i, particle_set[l], T, A_block, rho, a1, a2, 0.05);
      delete island_candidate;
      island_candidate = new Partition(particle_set[l]);
      best_split(isl_i, island_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda, false);
      island_obj = w[l]*total_log_post(island_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda*Entropy(l, island_candidate, particle_set,w);
      island_flag = Partition_Equal(island_candidate, particle_set[l]);
      if(island_flag == 0) try_local = false;
      if(island_obj > accepted_obj) accepted_obj = island_obj;
      
      if(try_local == true){
        get_local(local_si, particle_set[l], T, A_block, rho, a1, a2);
        delete local_candidate;
        local_candidate = new Partition(particle_set[l]);
        best_split(local_si, local_candidate, l, particle_set, w, ybar, total_ss, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda, false);
        local_obj = w[l] * total_log_post(local_candidate, total_ss, T, nu_sigma, lambda_sigma) + lambda * Entropy(l, local_candidate, particle_set, w);
        local_flag = Partition_Equal(local_candidate, particle_set[l]);
        if(local_obj > accepted_obj) accepted_obj = island_obj;
      }
      
      
      if(accepted_obj == spec_split_obj){
        // spectral split accepted
        if(spec_split_flag == 0) particle_set[l]->Copy_Partition(spec_split_candidate);
        else conv_counter++;
      } else if(accepted_obj == tail_split_obj){
        if(tail_split_flag == 0) particle_set[l]->Copy_Partition(tail_split_candidate);
        else conv_counter++;
      } else if(accepted_obj == km_split_obj){
        if(km_split_flag == 0) particle_set[l]->Copy_Partition(km_split_candidate);
        else conv_counter++;
      } else if(accepted_obj == merge_obj){
        if(merge_flag == 0) particle_set[l]->Copy_Partition(merge_candidate);
        else conv_counter++;
      } else if(accepted_obj == border_obj){
        if(border_flag == 0) particle_set[l]->Copy_Partition(border_candidate);
        else conv_counter++;
      } else if(accepted_obj == island_obj){
        if(island_flag == 0) particle_set[l]->Copy_Partition(island_candidate);
        else conv_counter++;
      } else if( (try_local == true) && (accepted_obj == local_obj)){
        if(local_flag == 0) particle_set[l]->Copy_Partition(local_candidate);
        else conv_counter++;
      }
    } // closes loop over the particle set
      // update the importance weights now
  
    update_w(particle_set, w, L, total_ss, T, nu_sigma, lambda_sigma, lambda);
    
    if(conv_counter == L){
      Rcpp::Rcout << "  None of the particles moved on this sweep!" << endl;
      flag = 1;
    }
    // compute the objective
    objective = lambda * Entropy(0, particle_set[0], particle_set, w);
    for(int l = 0; l < L; l++){
      objective += w[l] * total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
      tmp_log_like[l] = total_log_like(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
      tmp_log_prior[l] = total_log_prior(particle_set[l]);
      tmp_log_post[l] = total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
      for(int i = 0; i < n; i++) tmp_alpha_hat[l][i] = particle_set[l]->alpha_hat[i];
    }
    /*
     if(abs((objective - old_objective)/old_objective) < 0.01 * eps){
     Rcpp::Rcout << "[ensm_cluster_mean]: Objective has not increased much" << endl;
     //flag = 1;
     }
     */
    // update the trajectories
    objective_trajectory.push_back(objective);
    format_particle_set(particle_set, tmp_list);
    particle_set_trajectory.push_back(tmp_list);
    log_like_trajectory.push_back(tmp_log_like);
    log_prior_trajectory.push_back(tmp_log_prior);
    log_post_trajectory.push_back(tmp_log_post);
    alpha_trajectory.push_back(tmp_alpha_hat);
    
    Rcpp::Rcout << "   Number of stationary particles = " << conv_counter << endl;
    Rcpp::Rcout << "   objective = " << objective << "   old_objective = " << old_objective << "  %diff = " << 100.0 * abs( (objective - old_objective)/objective) << endl;
    
    iter++;
  } // closes the main loop
  int time2 = time(&tp);
  
  // Find the unique particles
  //std::vector<LPPartition> unik_particles;
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
  //Rcpp::List init_unik_particles_out;
  Rcpp::List unik_particles_out;
  //format_particle_set(init_unik_particles, init_unik_particles_out);
  format_particle_set(unik_particles, unik_particles_out);
  
  // Prepare the trajectory of the particle set for output
  // Also compute the log-likelihood, log-prior, and log-posterior for every particle visited by the system
  Rcpp::List particle_trajectory_out(particle_set_trajectory.size());
  arma::mat log_like_trajectory_out(log_like_trajectory.size(),L);
  arma::mat log_prior_trajectory_out(log_prior_trajectory.size(),L);
  arma::mat log_post_trajectory_out(log_post_trajectory.size(),L);
  arma::cube alpha_trajectory_out(n, L, alpha_trajectory.size());
 
  for(int ix = 0; ix < particle_set_trajectory.size(); ix++){
    particle_trajectory_out[ix] = particle_set_trajectory[ix];
    for(int l = 0; l < L; l++){
      log_like_trajectory_out(ix,l) = log_like_trajectory[ix][l];
      log_prior_trajectory_out(ix,l) = log_prior_trajectory[ix][l];
      log_post_trajectory_out(ix,l) = log_post_trajectory[ix][l];
      for(int i = 0; i < n; i++){
        alpha_trajectory_out(i, l, ix) = alpha_trajectory[ix][l][i];
      }
    }
  }

  /*
  Rcpp::Rcout << "Looking now at log-post trajectory" << endl;
  for(int ix = 0; ix < log_post_trajectory.size(); ix++){
    Rcpp::Rcout << "iteration " << ix << " :" ;
    for(int l = 0; l < L; l++) Rcpp::Rcout << " " << log_post_trajectory[ix][l];
    Rcpp::Rcout << endl;
  }
  */
  

  arma::mat alpha_hat_particle = arma::zeros<mat>(n, L_star); // estimates for each unique particle
  for(int l = 0; l < L_star; l++){
    for(int i = 0; i < n; i++){
      alpha_hat_particle(i,l) = unik_particles[l]->alpha_hat[i];
    }
  }
  arma::vec tmp_alpha = alpha_hat_particle.col(0);
  arma::mat alpha_hat_wtd = arma::zeros<mat>(n,L_star);
  double w_sum = 0.0;
  for(int l = 0; l < L_star; l++){
    arma::vec tmp_alpha = arma::zeros<vec>(n);
    //w_sum = 0.0;
    w_sum = 0.0;
    for(int ll = 0; ll <= l; ll++){
      tmp_alpha += pstar[ll] * alpha_hat_particle.col(ll);
      w_sum += pstar[ll];
    }
    //Rcpp::Rcout << "l = " << l << "w_sum = " << w_sum << endl;
    tmp_alpha /= w_sum;
    alpha_hat_wtd.col(l) = tmp_alpha;
  }
  
  Rcpp::List results;
  //results["init_unik_particles"] = init_unik_particles_out;
  //results["init_pstar"] = init_pstar;
  //results["init_counts"] = init_counts;
  //results["init_log_like"] = init_log_like;
  //results["init_log_prior"] = init_log_prior;
  //results["init_log_post"] = init_log_post;
  results["particles"] = unik_particles_out;
  results["pstar"] = pstar;
  results["counts"] = counts;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  results["alpha"] = alpha_hat_particle;
  results["time"] = time2 - time1;
  results["objective_trajectory"] = objective_trajectory;
  results["particle_trajectory"] = particle_trajectory_out;
  results["log_prior_trajectory"] = log_prior_trajectory_out;
  results["log_like_trajectory"] = log_like_trajectory_out;
  results["log_post_trajectory"] = log_post_trajectory_out;
  results["alpha_trajectory"] = alpha_trajectory_out;
  //results["alpha_hat_particle"] = alpha_hat_particle;
  //results["alpha_hat_wtd"] = alpha_hat_wtd;
  return(results);
}

