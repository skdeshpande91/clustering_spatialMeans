//
//  ensm_cluster_mean2.cpp
//
//  This only considers spectral splits, island moves, borders, and merges
//  Basically eliminates tail_splits and km_splits
//  Created by Sameer Deshpande on 13 Feb 019.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include <vector>

using namespace arma;
using namespace std;


using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List ensm_cluster_mean2(arma::vec ybar, const int T,  const arma::mat A_block, const int L, Rcpp::List gamma_init, const double a1 = 1.0, const double a2 = 1.0,
                             const double nu_sigma = 1, const double lambda_sigma = 1, const double rho = 0.99, const double lambda = 1.0, const double eta = 1.0,
                             const int max_iter = 10, const double eps = 1e-3, const double split_frac = 0.1)
{
  int n = ybar.size();
  Rcpp::Rcout << "n = " << n << endl;
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
  

  // for the main loop we need the following quantitites
  LPPartition spec_split_candidate = new Partition(gamma_0); // for the spectral splits
  LPPartition tail_split_candidate = new Partition(gamma_0); // for the tail splits
  LPPartition merge_candidate = new Partition(gamma_0); // for the merge candidates
  LPPartition border_candidate = new Partition(gamma_0); // for the border moves
  LPPartition island_candidate = new Partition(gamma_0); // for the island candidates
  
  split_info spec_si; // holds the information for spectral splits
  split_info tail_si; // holds the information for the tail splits
  split_info bi ; // holds information for border moves
  split_info isl_i; // holds information for island moves
  merge_info mi; // holds informatino for the merge moves
  
  double spec_split_obj = 0.0;
  double tail_split_obj = 0.0;
  double merge_obj = 0.0;
  double border_obj = 0.0;
  double island_obj = 0.0;
  double accepted_obj = 0.0;
  
  int spec_split_flag = 1; // essentially will hold the value of Partition_Equal(spec_split_candidate, particle_set[l])
  int tail_split_flag = 1; // holds value of Partition_Equal(tail_split_candidate, particle_set[l])
  int merge_flag = 1; // holds value of Partition_Equal(merge_candidate, particle_set[l]);
  int border_flag = 1; // holds value of Partition_Equal(border_candidate, particle_set[l]);
  int island_flag = 1; // holds value of Partition_Equal(island_flag, particle_set[l])


  int iter = 0;
  int conv_counter = 0; // count the number of particles that remain the same
  int flag = 0; //
  
  double old_objective = 0.0;
  double objective = 0.0;
  
  objective = Entropy(0, particle_set[0], particle_set, w);
  for(int l = 0; l < L; l++){
    objective += w[l] * total_log_post(particle_set[l], nu_sigma, lambda_sigma);
  }


  while((iter < max_iter) & (flag == 0)){
    Rcpp::Rcout << "[ensm_cluster_mean]: Starting iter " << iter << endl;
    old_objective = objective;
    
    conv_counter = 0; // counts the number of particles unchanged in our sweep
    for(int l = 0; l < L; l++){
      //Rcpp::Rcout << "Starting to update particle " << l << endl;
      spec_split_flag = 1;
      tail_split_flag = 1;
      border_flag = 1;
      merge_flag = 1;
      island_flag = 1;
      
      // Spectral splits
      //Rcpp::Rcout << "    Starting spectral split" << endl;
      get_spectral_split(spec_si, particle_set[l], T, A_block, rho, a1, a2, split_frac);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << spec_si.num_splits << " spectral splits" << endl;
      delete spec_split_candidate;
      spec_split_candidate = new Partition(particle_set[l]);
      best_split(spec_si, spec_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      spec_split_obj = w[l]*total_log_post(spec_split_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, spec_split_candidate, particle_set, w);
      
      spec_split_flag = Partition_Equal(spec_split_candidate, particle_set[l]);
      //Rcpp::Rcout << "  orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "    spec_split_obj = " << spec_split_obj << "    spec_split_flag =  " << spec_split_flag << endl;
      accepted_obj = spec_split_obj;
      
      // tail split
      get_tail_split(tail_si, particle_set[l], T, A_block, rho, a1, a2, split_frac);
      delete tail_split_candidate;
      tail_split_candidate = new Partition(particle_set[l]);
      best_split(tail_si, tail_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      tail_split_obj = w[l]*total_log_post(tail_split_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, tail_split_candidate, particle_set,w);
      tail_split_flag = Partition_Equal(tail_split_candidate, particle_set[l]);
      if(tail_split_obj > accepted_obj){
        accepted_obj = tail_split_obj;
      }
      // merges
      //Rcpp::Rcout << "    Starting merge" << endl;
      get_merge(mi, particle_set[l], A_block);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << mi.num_merges << " merge proposals" << endl;
      delete merge_candidate;
      merge_candidate = new Partition(particle_set[l]);
      best_merge(mi, merge_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      merge_obj = w[l]*total_log_post(merge_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, merge_candidate, particle_set, w);
      
      merge_flag = Partition_Equal(merge_candidate, particle_set[l]);
      //Rcpp::Rcout << "      orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "      merge_obj = " << merge_obj << "  merge_flag = " << merge_flag << endl;
      if(merge_obj > accepted_obj){
        accepted_obj = merge_obj;
      }
      // border
      //Rcpp::Rcout << "    Starting border" << endl;
      get_border(bi, particle_set[l], A_block);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << bi.num_splits << " border proposals" << endl;
      delete border_candidate;
      border_candidate = new Partition(particle_set[l]);
      //best_split(bi, border_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, a_sigma, nu_sigma, eta, lambda);
      //border_obj = w[l]*total_log_post(border_candidate, a_sigma, nu_sigma) + lambda*Entropy(l, border_candidate, particle_set,w);
      
      best_split(bi, border_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      border_obj = w[l]*total_log_post(border_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, border_candidate, particle_set,w);
      
      border_flag = Partition_Equal(border_candidate, particle_set[l]);
      //Rcpp::Rcout << "      orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "      border_obj = " << border_obj << "   border_flag = " << border_flag << endl;
      if(border_obj > accepted_obj){
        accepted_obj = border_obj;
      }
      //Rcpp::Rcout << "best border is " << endl;
      //border_candidate->Print_Partition();
      
      // island
      //Rcpp::Rcout << "    Starting island" << endl;
      get_island(isl_i, particle_set[l], A_block);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << isl_i.num_splits << " islands" << endl;
      delete island_candidate;
      island_candidate = new Partition(particle_set[l]);
      best_split(isl_i, island_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      island_obj = w[l]*total_log_post(island_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, island_candidate, particle_set,w);
      island_flag = Partition_Equal(island_candidate, particle_set[l]);
      
      //Rcpp::Rcout << "      orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "      island_obj = " << island_obj << "  island_flag = " << island_flag <<endl;

      if(island_obj > accepted_obj){
        accepted_obj = island_obj;
      }
      //Rcpp::Rcout << "  orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "  spec_split_obj = " << spec_split_obj << endl;
      //Rcpp::Rcout << "  merge_obj = " << merge_obj << endl;
      //Rcpp::Rcout << "  border_obj = " << border_obj << endl;
      //Rcpp::Rcout << "  island_obj = " << island_obj << endl;
      //Rcpp::Rcout << "  accepted_obj = " << accepted_obj << endl;

      if(accepted_obj == spec_split_obj){
        // spectral split accepted
        if(spec_split_flag == 0){
          particle_set[l]->Copy_Partition(spec_split_candidate);
        } else{
          // particle has not moved
          conv_counter++;
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved" << endl;
        }
      } else if(accepted_obj == tail_split_obj){
        if(tail_split_flag == 0){
          particle_set[l]->Copy_Partition(tail_split_candidate);
        } else{
          conv_counter++;
        }
      } else if(accepted_obj == merge_obj){
        if(merge_flag == 0){
          particle_set[l]->Copy_Partition(merge_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == merge_obj & merge_flag == 1" << endl;
          conv_counter++;
        }
      } else if(accepted_obj == border_obj){
        if(border_flag == 0){
          particle_set[l]->Copy_Partition(border_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == border_obj & border_flag == 1" << endl;
          conv_counter++;
        }
      } else if(accepted_obj == island_obj){
        if(island_flag == 0){
          particle_set[l]->Copy_Partition(island_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted obj == island_obj & island_flag == 1" << endl;
          conv_counter++;
        }
      }
      //Rcpp::Rcout << "Finished updating particle " << l << endl;
      //particle_set[l]->Print_Partition(nu_sigma, lambda_sigma);
      //Rcpp::Rcout << "log_post = " << total_log_post(particle_set[l]) << endl;
    } // closes loop over the particle set
    // update the importance weights now
    //Rcpp::Rcout << "About to update w" << endl;
    //update_w(particle_set, w, L, a_sigma, nu_sigma,lambda);
    update_w(particle_set, w, L, nu_sigma, lambda_sigma, lambda);
    //for(int l = 0;l < L; l++){
    //  Rcpp::Rcout << "w[" << l << "] = " << w[l] << endl;
    //}
    
    if(conv_counter == L){
      Rcpp::Rcout << "None of the particles moved on this sweep!" << endl;
      flag = 1;
    }
    // compute the objective
    objective = Entropy(0, particle_set[0], particle_set, w);
    for(int l = 0; l < L; l++){
      //objective += w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma);
      objective += w[l] * total_log_post(particle_set[l], nu_sigma, lambda_sigma);

    }
    if(abs((objective - old_objective)/old_objective) < 0.01 * eps){
      Rcpp::Rcout << "[ensm_cluster_mean]: Objective has not increased much" << endl;
      flag = 1;
    }
    Rcpp::Rcout << "   conv_counter = " << conv_counter << endl;
    Rcpp::Rcout << "   objective = " << objective << "   old_objective = " << old_objective << "  %diff = " << abs( (objective - old_objective)/objective) << endl;
    
    iter++;
    //Rcpp::Rcout << "Particle 0 is now" << endl;
    //particle_set[0]->Print_Partition(a_sigma, nu_sigma);
  }


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
    log_like[l] = total_log_like(unik_particles[l], nu_sigma, lambda_sigma);
    log_prior[l] = total_log_prior(unik_particles[l]);
    log_post[l] = total_log_post(unik_particles[l], nu_sigma, lambda_sigma);

  }
  
  arma::vec tmp_vec = arma::zeros<vec>(1);
  Rcpp::List tmp_list;
  Rcpp::List init_unik_particles_out(init_L_star);
  Rcpp::List unik_particles_out(L_star);
  Rcpp::NumericVector output_vec;
  
  for(int l = 0; l < init_L_star; l++){
    tmp_list = List(init_unik_particles[l]->K);
    for(int k = 0; k < init_unik_particles[l]->K; k++){
      tmp_vec = arma::zeros<vec>(init_unik_particles[l]->cluster_config[k]);
      for(int i = 0; i < init_unik_particles[l]->cluster_config[k]; i++){
        tmp_vec(i) = init_unik_particles[l]->clusters[k][i] + 1; // remember R is 1-indexed
      }
      output_vec = Rcpp::wrap(arma::sort(tmp_vec, "ascend"));
      output_vec.attr("dim") = R_NilValue;
      tmp_list[k] = output_vec;
    }
    init_unik_particles_out[l] = tmp_list;
  }
  
  for(int l = 0; l < L_star; l++){
    tmp_list = List(unik_particles[l]->K);
    for(int k = 0; k < unik_particles[l]->K; k++){
      tmp_vec = arma::zeros<vec>(unik_particles[l]->cluster_config[k]);
      for(int i = 0; i < unik_particles[l]->cluster_config[k];i++){
        tmp_vec(i) = unik_particles[l]->clusters[k][i] + 1; // remember R is 1-indexed
      }
      output_vec = Rcpp::wrap(arma::sort(tmp_vec, "ascend"));
      output_vec.attr("dim") = R_NilValue;
      tmp_list[k] = output_vec;
    }
    unik_particles_out[l] = tmp_list;
  }
  
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
  results["init_unik_particles"] = init_unik_particles_out;
  results["init_pstar"] = init_pstar;
  results["init_counts"] = init_counts;
  results["init_log_like"] = init_log_like;
  results["init_log_prior"] = init_log_prior;
  results["init_log_post"] = init_log_post;
  results["particles"] = unik_particles_out;
  results["pstar"] = pstar;
  results["counts"] = counts;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  results["alpha_hat_particle"] = alpha_hat_particle;
  results["alpha_hat_wtd"] = alpha_hat_wtd;
  return(results);
}
