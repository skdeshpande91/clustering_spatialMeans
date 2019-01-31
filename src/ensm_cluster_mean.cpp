//
//  ensm_cluster_mean.cpp
//  
//
//  Created by Sameer Deshpande on 10/29/18.
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
Rcpp::List ensm_cluster_mean(arma::vec ybar, const int T,  const arma::mat A_block, const int L, Rcpp::List gamma_init, const double a1 = 1.0, const double a2 = 1.0,
                             const double nu_sigma = 1, const double lambda_sigma = 1, const double rho = 0.99, const double lambda = 1.0, const double eta = 1.0,
                             const int max_iter = 10, const double eps = 1e-3, const double split_frac = 0.1)
{
  int n = ybar.size();
  Rcpp::Rcout << "n = " << n << endl;
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::Rcout << "Created gamma_0" << endl;
  
  
  // We should start with everything in one partition. then we can run split operations (e.g. k-means or spectral clustering or tail_splits)
  // Then we initialize particle set with at least two copies of each thing found
  // Then we sample the rest
  
  
  
  //for(int k = 0; k < gamma_0->K; k++){
  //  Rcpp::Rcout << gamma_0->log_det_Omegay[k] << " ";
  //}
  //Rcpp::Rcout << endl;
  

  // create a particle set
  std::vector<LPPartition> particle_set(L);
  std::vector<double> w(L);
  for(int l = 0; l < L; l++){
    particle_set[l] = new Partition(gamma_0);
    w[l] = 1.0/( (double) L);
  }

  // [SKD]: 12 November 2018 -- at this point, run spectral clustering starting from cluster with everything togeter.
  //        then sort and do importance sampling. Copy the code from BART stuff for random number generation if necessary
  

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
    //init_log_like[l] = total_log_like(init_unik_particles[l], a_sigma, nu_sigma);
    init_log_like[l] = total_log_like(init_unik_particles[l], nu_sigma, lambda_sigma);
    init_log_prior[l] = total_log_prior(init_unik_particles[l]);
    //init_log_post[l] = total_log_post(init_unik_particles[l], a_sigma, nu_sigma);
    init_log_post[l] = total_log_post(init_unik_particles[l], nu_sigma, lambda_sigma);

  }
  

  // for the main loop we need the following quantitites
  LPPartition spec_split_candidate = new Partition(gamma_0); // for the spectral splits
  LPPartition tail_split_candidate = new Partition(gamma_0); // for the tail splits
  LPPartition km_split_candidate = new Partition(gamma_0); // for the K-means splits
  LPPartition merge_candidate = new Partition(gamma_0); // for the merge candidates
  LPPartition border_candidate = new Partition(gamma_0); // for the border moves
  LPPartition island_candidate = new Partition(gamma_0); // for the island candidates
  
  //LPPartition accepted_candidate = new Partition(gamma_0); // the actual canidate that is accepte
  
  split_info spec_si; // holds the information for spectral splits
  split_info tail_si; // holds information for the tail splits
  split_info km_si; // holds info for the k-mean splits
  split_info bi ; // holds information for border moves
  split_info isl_i; // holds information for island moves
  //split_info mi; // holds information for merge moves
  merge_info mi; // holds informatino for the merge moves
  
  double spec_split_obj = 0.0;
  double tail_split_obj = 0.0;
  double km_split_obj = 0.0;
  double merge_obj = 0.0;
  double border_obj = 0.0;
  double island_obj = 0.0;
  double accepted_obj = 0.0;
  
  int spec_split_flag = 1; // essentially will hold the value of Partition_Equal(spec_split_candidate, particle_set[l])
  int tail_split_flag = 1; // holds value of Partition_Equal(split_candidate, particle_set[l])
  int km_split_flag = 1; // holds value of Partition_Equal(km_split_candidate, particle_set[l]);
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
    //objective += w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma);
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
      km_split_flag = 1;
      border_flag = 1;
      merge_flag = 1;
      island_flag = 1;
      
      // Spectral splits
      //Rcpp::Rcout << "    Starting spectral split" << endl;
      get_spectral_split(spec_si, particle_set[l], T, A_block, rho, a1, a2, split_frac);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << spec_si.num_splits << " spectral splits" << endl;
      delete spec_split_candidate;
      spec_split_candidate = new Partition(particle_set[l]);
      //best_split(spec_si, spec_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, a_sigma, nu_sigma, eta, lambda);
      //spec_split_obj = w[l]*total_log_post(spec_split_candidate, a_sigma, nu_sigma) + lambda * Entropy(l, spec_split_candidate, particle_set, w);
      
      best_split(spec_si, spec_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma, eta, lambda);
      spec_split_obj = w[l]*total_log_post(spec_split_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, spec_split_candidate, particle_set, w);
      
      spec_split_flag = Partition_Equal(spec_split_candidate, particle_set[l]);
      //Rcpp::Rcout << "  orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "    spec_split_obj = " << spec_split_obj << "    spec_split_flag =  " << spec_split_flag << endl;
      accepted_obj = spec_split_obj;
      
      // tail splits
      //Rcpp::Rcout << "    Starting tail split" << endl;
      get_tail_split(tail_si, particle_set[l], T, A_block, rho, a1, a2, split_frac);
      //Rcpp::Rcout <<"[ensm_cluster_mean]: Got " << tail_si.num_splits << " tail splits" << endl;
      delete tail_split_candidate;
      tail_split_candidate = new Partition(particle_set[l]);
      //best_split(tail_si, tail_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2,a_sigma, nu_sigma, eta, lambda);
      //tail_split_obj = w[l]*total_log_post(tail_split_candidate, a_sigma, nu_sigma) + lambda*Entropy(l, tail_split_candidate, particle_set, w);
      
      best_split(tail_si, tail_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2,nu_sigma, lambda_sigma, eta, lambda);
      tail_split_obj = w[l]*total_log_post(tail_split_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, tail_split_candidate, particle_set, w);
      
      tail_split_flag = Partition_Equal(tail_split_candidate, particle_set[l]);
      //Rcpp::Rcout << "      orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "      tail_split_obj = " << tail_split_obj << "   tail_split_flag = " << tail_split_flag << endl;

      if(tail_split_obj > accepted_obj){
        accepted_obj = tail_split_obj;
      }
      // K-means splits
      //Rcpp::Rcout << "    Starting KM split" << endl;
      get_km_split(km_si, particle_set[l], T, A_block, rho, a1, a2, 5*split_frac);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << km_si.num_splits << " k-means splits" << endl;
      delete km_split_candidate;
      km_split_candidate = new Partition(particle_set[l]);
      //best_split(km_si, km_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, a_sigma, nu_sigma,eta, lambda);
      //km_split_obj = w[l]*total_log_post(km_split_candidate, a_sigma, nu_sigma) + lambda*Entropy(l, km_split_candidate, particle_set, w);
      
      best_split(km_si, km_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, nu_sigma, lambda_sigma,eta, lambda);
      km_split_obj = w[l]*total_log_post(km_split_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, km_split_candidate, particle_set, w);
      
      km_split_flag = Partition_Equal(km_split_candidate, particle_set[l]);
      //Rcpp::Rcout << "      orig_obj = " << w[l] * total_log_post(particle_set[l], a_sigma, nu_sigma) + lambda*Entropy(l, particle_set[l], particle_set, w) << endl;
      //Rcpp::Rcout << "      km_split_obj = " << km_split_obj << "    km_split_flag = " << km_split_flag << endl;
      if(km_split_obj > accepted_obj){
        accepted_obj = km_split_obj;
      }
      // merges
      //Rcpp::Rcout << "    Starting merge" << endl;
      get_merge(mi, particle_set[l], A_block);
      //Rcpp::Rcout << "[ensm_cluster_mean]: Got " << mi.num_merges << " merge proposals" << endl;
      delete merge_candidate;
      merge_candidate = new Partition(particle_set[l]);
      //best_split(mi, merge_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, a_sigma, nu_sigma, eta, lambda);
      //best_merge(mi, merge_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, a_sigma, nu_sigma, eta, lambda);
      //merge_obj = w[l]*total_log_post(merge_candidate, a_sigma, nu_sigma) + lambda*Entropy(l, merge_candidate, particle_set, w);
      
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
      //best_split(isl_i, island_candidate, l, particle_set, w, ybar, T, A_block, rho, a1, a2, a_sigma, nu_sigma, eta, lambda);
      //island_obj = w[l]*total_log_post(island_candidate, a_sigma, nu_sigma) + lambda*Entropy(l, island_candidate, particle_set,w);
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
      //Rcpp::Rcout << "  tail_split_obj = " << tail_split_obj << endl;
      //Rcpp::Rcout << "  km_split_obj = " << km_split_obj << endl;
      //Rcpp::Rcout << "  merge_obj = " << merge_obj << endl;
      //Rcpp::Rcout << "  border_obj = " << border_obj << endl;
      //Rcpp::Rcout << "  island_obj = " << island_obj << endl;
      //Rcpp::Rcout << "  accepted_obj = " << accepted_obj << endl;

      if(accepted_obj == spec_split_obj){
        // spectral split accepted
        if(spec_split_flag == 0){
          //delete particle_set[l];
          //particle_set[l] = new Partition(spec_split_candidate);
          //Rcpp::Rcout << "  spectral split accepted" << endl;
          particle_set[l]->Copy_Partition(spec_split_candidate);
        } else{
          // particle has not moved
          conv_counter++;
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved" << endl;
        }
      } else if(accepted_obj == tail_split_obj){
        if(tail_split_flag == 0){
          //delete particle_set[l];
          //particle_set[l] = new Partition(tail_split_candidate);
          //Rcpp::Rcout << "  Tail split accepted" << endl;
          particle_set[l]->Copy_Partition(tail_split_candidate);
        } else{
          // control should never reach here. If particle doesn't move, then all of the *_obj will be equal. We only let accepted_obj != spec_split_obj if some other proposal yields higher objective. But since we always consider particle_set[l] among all of the proposed moves, we can't update accepted_obj unless the update is non-trivial
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved... but control reached bad spot with accepted_obj == tail_split_obj & tail_split_flag == 1" << endl;
          conv_counter++;
        }
      } else if(accepted_obj == km_split_obj){
        if(km_split_flag == 0){
          //delete particle_set[l];
          //particle_set[l] = new Partition(km_split_candidate);
          //Rcpp::Rcout << "  KM split accepted" << endl;
          particle_set[l]->Copy_Partition(km_split_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == km_split_obj & km_split_flag == 1" << endl;
          // particle has not moved but control should never reach here
          conv_counter++;
        }
      } else if(accepted_obj == merge_obj){
        if(merge_flag == 0){
          //delete particle_set[l];
          //particle_set[l] = new Partition(merge_candidate);
          //Rcpp::Rcout << "  merge split accepted" << endl;
          particle_set[l]->Copy_Partition(merge_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == merge_obj & merge_flag == 1" << endl;
          conv_counter++;
        }
      } else if(accepted_obj == border_obj){
        if(border_flag == 0){
          //delete particle_set[l];
          //particle_set[l] = new Partition(border_candidate);
          //Rcpp::Rcout << "  border move accepted" << endl;
          particle_set[l]->Copy_Partition(border_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == border_obj & border_flag == 1" << endl;
          conv_counter++;
        }
      } else if(accepted_obj == island_obj){
        if(island_flag == 0){
          //delete particle_set[l];
          //particle_set[l] = new Partition(island_candidate);
          //Rcpp::Rcout << "  island move accepted" << endl;
          particle_set[l]->Copy_Partition(island_candidate);
        } else{
          //Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted obj == island_obj & island_flag == 1" << endl;
          conv_counter++;
        }
      }
      //Rcpp::Rcout << "Finished updating particle " << l << endl;
      //particle_set[l]->Print_Partition(a_sigma, nu_sigma);
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
    //log_like[l] = total_log_like(unik_particles[l], a_sigma, nu_sigma);
    log_like[l] = total_log_like(unik_particles[l], nu_sigma, lambda_sigma);
    log_prior[l] = total_log_prior(unik_particles[l]);
    //log_post[l] = total_log_post(unik_particles[l], a_sigma, nu_sigma);
    log_post[l] = total_log_post(unik_particles[l], nu_sigma, lambda_sigma);
    
    
    // 14 January Update: we are getting some really strange results with the log-likelihood being computed at the end here
    //Rcpp::Rcout << "unik partition" << l << endl;
    //unik_particles[l]->Print_Partition(a_sigma, nu_sigma);
    //Rcpp::Rcout << "log_like[" << l << "] = " << log_like[l] << endl;
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

  
  /*
  Rcpp::Rcout << "Finished with " << particle_map.size() << " unique particles" << endl;
  int L_star = unik_particles.size();
  arma::vec log_like = arma::zeros<vec>(L_star);
  arma::vec log_prior = arma::zeros<vec>(L_star);
  arma::vec log_post = arma::zeros<vec>(L_star);



  arma::vec tmp_vec = arma::zeros<vec>(1);
  Rcpp::List tmp_list; // for the final particle set
  Rcpp::List unik_particles_out(L_star);
  Rcpp::NumericVector output_vec;
  
  Rcpp::Rcout << "Preparing output for final particle set" << endl;
  for(int l = 0; l < L_star; l++){
    log_like[l] = total_log_like(unik_particles[l], a_sigma, nu_sigma);
    log_prior[l] = total_log_prior(unik_particles[l]);
    log_post[l] = total_log_post(unik_particles[l], a_sigma, nu_sigma);
    
    tmp_list = List(unik_particles[l]->K);
    for(int k = 0; k < unik_particles[l]->K; k++){
      tmp_vec = arma::zeros<vec>(unik_particles[l]->cluster_config[k]);
      for(int i = 0; i < unik_particles[l]->cluster_config[k]; i++){
        tmp_vec[i] = unik_particles[l]->clusters[k][i]+1; // remember that in R we are 1-indexed
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
    w_sum = 0.0;
    for(int ll = 0; ll < l; ll++){
      tmp_alpha += p_star[ll] * alpha_hat_particle.col(ll);
      w_sum += p_star[ll];
    }
    tmp_alpha /= w_sum;
    alpha_hat_wtd.col(l) = tmp_alpha;
  }

  // clean up more memory
  delete gamma_0;
  //for(int l = 0; l > init_L_star; l++){
  //  delete init_particle_set[l];
  //}
  //init_particle_set.clear();
  for(int l = 0; l < L; l++){
    delete particle_set[l];
  }
  particle_set.clear();

  for(int l = 0; l < L_star; l++){
    delete unik_particles[l];
  }
  unik_particles.clear();
  
  delete spec_split_candidate;
  delete tail_split_candidate;
  delete km_split_candidate;
  delete merge_candidate;
  delete border_candidate;
  delete island_candidate;
  

  Rcpp::List results;
  //results["ybar"] = ybar;
  //results["init_particle_set"] = init_particles_out;
  //results["init_w"] = init_w;
  //results["init_counts"] = init_counts;
  //results["init_log_like"] = init_log_like;
  //results["init_log_prior"] = init_log_prior;
  //results["init_log_post"] = init_log_post;
  results["particle_set"] = unik_particles_out;
  results["w"] = p_star;
  results["counts"] = counts;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  results["alpha_hat_particle"] = alpha_hat_particle;
  results["alpha_hat_wtd"] = alpha_hat_wtd;
   */


}
/*
 // Try a border move
 border_info bi;
 double border_frac = 0.01;
 Rcpp::Rcout << "created bi" << endl;
 get_border(bi, gamma_0, A_block, border_frac);
 Rcpp::Rcout << "Proposing " << bi.num_borders << "  border candidates:" << endl;
 
 LPPartition border_candidate = new Partition(gamma_0);
 best_border(border_candidate, 0, bi, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 Rcpp::Rcout << "Best border is:" << endl;
 border_candidate->Print_Partition();
 
 
 double split_frac = 0.05;
 split_info si;
 get_new_split(si, gamma_0, T, A_block, rho, a, split_frac);
 Rcpp::Rcout << "Proposed " << si.num_splits << " potential splits" << endl;
 LPPartition split_candidate = new Partition(gamma_0);
 best_split(si, split_candidate, 0, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 Rcpp::Rcout << "Best new split is " << endl;
 split_candidate->Print_Partition();
 
 split_info spec_si;
 get_spectral_split(spec_si, gamma_0, T, A_block, rho, a, split_frac);
 LPPartition spec_split_candidate = new Partition(gamma_0);
 best_split(spec_si, spec_split_candidate, 0, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 Rcpp::Rcout << "Best spectral Split is: " << endl;
 spec_split_candidate->Print_Partition();
 
 // try a merge move
 split_info mi;
 get_merge(mi, gamma_0, A_block);
 Rcpp::Rcout << "Proposed " << mi.num_splits << "merges" << endl;
 for(int ix = 0; ix < mi.num_splits; ix++){
 Rcpp::Rcout << "Merge cluster " << mi.split_k[ix] << " into " << mi.nearest_neighbor[ix][0] << endl;
 }
 
 LPPartition merge_candidate = new Partition(gamma_0);
 best_split(mi, merge_candidate, 0, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 Rcpp::Rcout << "best merge candidate is : " << endl;
 merge_candidate->Print_Partition();

 
 split_info local_i;
 get_local(local_i, gamma_0, A_block);
 Rcpp::Rcout << "Proposed " << local_i.num_splits << " local candidates" << endl;
 // print out the first local candidate
 for(int new_k = 0; new_k < local_i.new_clusters[0].size(); new_k++){
 Rcpp::Rcout << "  new cluster " << new_k << " of size " << local_i.new_clusters[0][new_k].size() << " and neighbor " << local_i.nearest_neighbor[0][new_k] << " : " << endl;
 for(int i = 0; i < local_i.new_clusters[0][new_k].size(); i++){
 Rcpp::Rcout << local_i.new_clusters[0][new_k][i] << " ";
 }
 Rcpp::Rcout << endl;
 }
 // [SKD]: get local is super slow.
 
 
 get_km_split(km_si, particle_set[l], T, A_block, rho, a, split_frac);
 Rcpp::Rcout << "[ensm_cluster_mean]: Got " << km_si.num_splits << " k-means splits" << endl;
 
 // print out the km splits
 for(int ix = 0; ix < km_si.num_splits; ix++){
 Rcpp::Rcout << "  Splitting cluster " << km_si.split_k[ix] <<  "  into " << km_si.new_clusters[ix].size() << "  new clusters:" << endl;
 for(int new_k = 0; new_k < km_si.new_clusters[ix].size(); new_k++){
 Rcpp::Rcout << "  new cluster " << new_k << " of size " << km_si.new_clusters[ix][new_k].size() << "  and neighbor " << km_si.nearest_neighbor[ix][new_k] << endl;
 for(int ii = 0; ii < km_si.new_clusters[ix][new_k].size(); ii++){
 Rcpp::Rcout << km_si.new_clusters[ix][new_k][ii] << " ";
 }
 Rcpp::Rcout << endl;
 }
 }
 */
/*
 
 // start with spectral splits
 int l = 0;
 
 
 get_spectral_split(spec_si, particle_set[l], T, A_block, rho, a, split_frac);
 Rcpp::Rcout << "[ensm_cluster_mean]: Got " << spec_si.num_splits << " spectral splits" << endl;
 delete spec_split_candidate;
 spec_split_candidate = new Partition(particle_set[l]);
 best_split(spec_si, spec_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 spec_split_obj = w[l]*total_log_post(spec_split_candidate) + lambda * Entropy(l, spec_split_candidate, particle_set, w);
 spec_split_flag = Partition_Equal(spec_split_candidate, particle_set[l]);
 
 accepted_obj = spec_split_obj;
 //Rcpp::Rcout << "best spectral split is " << endl;
 //spec_split_candidate->Print_Partition();
 
 
 
 // tail splits
 get_tail_split(tail_si, particle_set[l], T, A_block, rho, a, split_frac);
 Rcpp::Rcout <<"[ensm_cluster_mean]: Got " << tail_si.num_splits << " tail splits" << endl;
 delete tail_split_candidate;
 tail_split_candidate = new Partition(particle_set[l]);
 best_split(tail_si, tail_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 tail_split_obj = w[l]*total_log_post(tail_split_candidate) + lambda*Entropy(l, tail_split_candidate, particle_set, w);
 tail_split_flag = Partition_Equal(tail_split_candidate, particle_set[l]);
 if(tail_split_obj > accepted_obj){
 accepted_obj = tail_split_obj;
 }
 //Rcpp::Rcout << "best tail split is " << endl;
 //tail_split_candidate->Print_Partition();
 
 // K-means splits
 get_km_split(km_si, particle_set[l], T, A_block, rho, a, 5*split_frac);
 Rcpp::Rcout << "[ensm_cluster_mean]: Got " << km_si.num_splits << " k-means splits" << endl;
 delete km_split_candidate;
 km_split_candidate = new Partition(particle_set[l]);
 best_split(km_si, km_split_candidate, l, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 km_split_obj = w[l]*total_log_post(km_split_candidate) + lambda*Entropy(l, km_split_candidate, particle_set, w);
 km_split_flag = Partition_Equal(km_split_candidate, particle_set[l]);
 if(km_split_obj > accepted_obj){
 accepted_obj = km_split_obj;
 }
 //Rcpp::Rcout << "best km split is " << endl;
 // km_split_candidate->Print_Partition();
 
 // merges
 get_merge(mi, particle_set[l], A_block);
 Rcpp::Rcout << "[ensm_cluster_mean]: Got " << mi.num_splits << " merge proposals" << endl;
 delete merge_candidate;
 merge_candidate = new Partition(particle_set[l]);
 best_split(mi, merge_candidate, l, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 merge_obj = w[l]*total_log_post(merge_candidate) + lambda*Entropy(l, merge_candidate, particle_set, w);
 merge_flag = Partition_Equal(merge_candidate, particle_set[l]);
 if(merge_obj > accepted_obj){
 accepted_obj = merge_obj;
 }
 //Rcpp::Rcout << "best merge is " << endl;
 //merge_candidate->Print_Partition();
 
 // border
 get_border(bi, particle_set[l], A_block);
 Rcpp::Rcout << "[ensm_cluster_mean]: Got " << bi.num_splits << " border proposals" << endl;
 delete border_candidate;
 border_candidate = new Partition(particle_set[l]);
 best_split(bi, border_candidate, l, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 border_obj = w[l]*total_log_post(border_candidate) + lambda*Entropy(l, border_candidate, particle_set,w);
 border_flag = Partition_Equal(border_candidate, particle_set[l]);
 if(border_obj > accepted_obj){
 accepted_obj = border_obj;
 }
 //Rcpp::Rcout << "best border is " << endl;
 //border_candidate->Print_Partition();
 
 // island
 get_island(isl_i, particle_set[l], A_block);
 Rcpp::Rcout << "[ensm_cluster_mean]: Got " << isl_i.num_splits << " islands" << endl;
 delete island_candidate;
 island_candidate = new Partition(particle_set[l]);
 best_split(isl_i, island_candidate, l, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 island_obj = w[l]*total_log_post(island_candidate) + lambda*Entropy(l, island_candidate, particle_set,w);
 island_flag = Partition_Equal(island_candidate, particle_set[l]);
 if(island_obj > accepted_obj){
 accepted_obj = island_obj;
 }
 //Rcpp::Rcout << "best island is " << endl;
 //island_candidate->Print_Partition();
 
 Rcpp::Rcout << "accepted_obj = " << accepted_obj << endl;
 Rcpp::Rcout << "spec_split_obj = " << spec_split_obj << endl;
 Rcpp::Rcout << "tail_split_obj = " << tail_split_obj << endl;
 Rcpp::Rcout << "km_split_obj = " << km_split_obj << endl;
 Rcpp::Rcout << "merge_obj = " << merge_obj << endl;
 Rcpp::Rcout << "border_obj = " << border_obj << endl;
 Rcpp::Rcout << "island_obj = " << island_obj << endl;
 
 
 // figure out which proposal was accepted
 if(accepted_obj == spec_split_obj){
 // spectral split accepted
 if(spec_split_flag == 0){
 delete particle_set[l];
 particle_set[l] = new Partition(spec_split_candidate);
 } else{
 // particle has not moved
 Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved" << endl;
 }
 } else if(accepted_obj == tail_split_obj){
 if(tail_split_flag == 0){
 delete particle_set[l];
 particle_set[l] = new Partition(tail_split_candidate);
 } else{
 // control should never reach here. If particle doesn't move, then all of the *_obj will be equal. We only let accepted_obj != spec_split_obj if some other proposal yields higher objective. But since we always consider particle_set[l] among all of the proposed moves, we can't update accepted_obj unless the update is non-trivial
 Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved... but control reached bad spot with accepted_obj == tail_split_obj & tail_split_flag == 1" << endl;
 }
 } else if(accepted_obj == km_split_obj){
 if(km_split_flag == 0){
 delete particle_set[l];
 particle_set[l] = new Partition(km_split_candidate);
 } else{
 Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == km_split_obj & km_split_flag == 1" << endl;
 }
 } else if(accepted_obj == merge_obj){
 if(merge_flag == 0){
 delete particle_set[l];
 particle_set[l] = new Partition(merge_candidate);
 } else{
 Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == merge_obj & merge_flag == 1" << endl;
 }
 } else if(accepted_obj == border_obj){
 if(border_flag == 0){
 delete particle_set[l];
 particle_set[l] = new Partition(border_candidate);
 } else{
 Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted_obj == border_obj & border_flag == 1" << endl;
 }
 } else if(accepted_obj == island_obj){
 if(island_flag == 0){
 delete particle_set[l];
 particle_set[l] = new Partition(island_candidate);
 } else{
 Rcpp::Rcout << "[ensm_cluster_mean]: particle has not moved ... but control reached bad spot with accepted obj == island_obj & island_flag == 1" << endl;
 }
 }
 
 particle_set[l]->Print_Partition();
 
*/



/*
 get_spectral_split(si, gamma_0, T, A_block, rho, a, split_frac);
 
 Rcpp::Rcout << "Proposed " << si.num_splits << " potential splits" << endl;
 for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
 Rcpp::Rcout << "   Split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " new clusters" << endl;
 for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
 Rcpp::Rcout << " new cluster " << nc_ix << "  whose nearest neighbor is  " << si.nearest_neighbor[split_ix][nc_ix] <<" : ";
 for(int i = 0; i < si.new_clusters[split_ix][nc_ix].size(); i++){
 Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][i] << " " ;
 }
 Rcpp::Rcout << endl;
 }
 }
 
 //
 Rcpp::Rcout << "Finding the best split now" << endl;
 LPPartition split_candidate = new Partition(gamma_0);
 best_split(si, split_candidate, 0, particle_set, w, ybar, T, A_block, rho, a, eta, lambda);
 split_candidate->Print_Partition();
 */

/*
 // we are going to try Split_Merge now starting with gamma_4
 std::vector<std::vector<int> > new_clusters(3);
 // will split cluster labelled 1.
 for(int j = 10; j < 13; j++){
 for(int i = 0; i < 10; i++){
 new_clusters[0].push_back(i + 20*j);
 }
 }
 for(int j = 13; j < 15; j++){
 for(int i = 0; i < 10; i++){
 new_clusters[1].push_back(i + 20*j);
 }
 }
 for(int j = 15; j < 20; j++){
 for(int i = 0; i < 10; i++){
 new_clusters[2].push_back(i + 20*j);
 }
 }
 
 std::vector<int> k_star;
 k_star.push_back(0);
 k_star.push_back(3);
 k_star.push_back(3);
 
 Rcpp::Rcout << "Ready to try to split" << endl;
 gamma_0->Split_Merge(1, new_clusters, k_star, ybar, T, A_block, rho, a, eta);
 gamma_0->Print_Partition();
 */
