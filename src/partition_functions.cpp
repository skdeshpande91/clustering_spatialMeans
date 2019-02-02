/*
 * partition_functions.cpp
 *
 *  Created on: Dec 29, 2017
 *      Author: Sameer
 */
//#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <algorithm>
#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"



// Compare two partitions and see if they are equal
int Partition_Equal(Partition *partition1, Partition *partition2){
  int flag = 1;
    // simpler to just compare pairwise allocation
  for(int i = 0 ; i < partition1->nObs;i++){
    for(int j = 0; j < partition1->nObs;j++){
      if(partition1->pairwise_assignment(i,j) != partition2->pairwise_assignment(i,j)){
        //Rcpp::Rcout << "[Partition_Equal]: i = " << i << " j = " << j << endl;
        flag = 0;
        break;
      }
    }
    if(flag == 0) break;
  }
  return flag;
}

void get_unik_particles(std::vector<std::vector<int> > &particle_map, std::vector<double> &p_star, std::vector<int> &counts, std::vector<LPPartition> particle_set, std::vector<double> w)
{
  particle_map.clear();
  p_star.clear();
  counts.clear();
  std::vector<std::vector<int> > tmp_particle_map(1, std::vector<int>(1,0)); // first particle always considered unique
  std::vector<double> tmp_w(1, w[0]);
  int counter = 0;
  for(int l = 1; l < particle_set.size(); l++){
    counter = 0;
    for(int ul = 0; ul < tmp_particle_map.size(); ul++){
      if(Partition_Equal(particle_set[l], particle_set[tmp_particle_map[ul][0]]) == 1){
        tmp_particle_map[ul].push_back(l);
        tmp_w[ul] += w[l];
      } else{
        counter++;
      }
    }
    if(counter == tmp_particle_map.size()){
      tmp_particle_map.push_back(std::vector<int>(1,l));
      tmp_w.push_back(w[l]);
    }
  }
  
  arma::vec w_vec = arma::zeros<vec>(tmp_particle_map.size());
  arma::uvec w_indices(tmp_particle_map.size());
  for(int l = 0; l < tmp_particle_map.size(); l++){
    w_vec(l) = tmp_w[l];
  }
  w_indices = arma::sort_index(w_vec, "descend");
  for(int l = 0; l < tmp_particle_map.size(); l++){
    particle_map.push_back(tmp_particle_map[w_indices(l)]);
    p_star.push_back(tmp_w[w_indices(l)]);
    counts.push_back(tmp_particle_map[w_indices(l)].size());
  }
  
}


/*
void get_unik_particles(std::vector<LPPartition> &unik_particles, std::vector<double> &p_star, std::vector<int> &counts, std::vector<LPPartition> particle_set, std::vector<double> w)
{
  unik_particles.clear();
  p_star.clear();
  counts.clear();
  LPPartition tmp_partition = new Partition(particle_set[0]);
  std::vector<LPPartition> tmp_unik_particles;
  tmp_unik_particles.push_back(tmp_partition);
  
  //std::vector<LPPartition> tmp_unik_particles(1,particle_set[0]);
  std::vector<double> tmp_p_star(1, w[0]);
  std::vector<int> tmp_counts(1,1);

  int num_unik_particles = 1;
  int counter = 0;
  int L = particle_set.size();
  for(int l = 1; l < L; l++){
    //Rcpp::Rcout << "[get_unik_particles]: l = " << l << endl;
    counter = 0;
    for(int ul = 0; ul < num_unik_particles; ul++){
      if(Partition_Equal(particle_set[l], tmp_unik_particles[ul]) == 1){
        //Rcpp::Rcout << "    particle " << l << " equals unik particle " << ul << endl;
        tmp_p_star[ul] += w[l];
        tmp_counts[ul] ++;
        break;
      } else{
        counter++;
      }
    }
    if(counter == num_unik_particles){
      //Rcpp::Rcout << "[get_unik_particles]: Found a unique particle!" << endl;
      //tmp_unik_particles.push_back(particle_set[l]);
      delete tmp_partition;
      tmp_partition = new Partition(particle_set[l]);
      tmp_unik_particles.push_back(tmp_partition);
      tmp_p_star.push_back(w[l]);
      tmp_counts.push_back(1);
      num_unik_particles++;
    }
  }
  //Rcpp::Rcout << "[get_unik_particles]: Got " << num_unik_particles << " unique particles" << endl;
  //Rcpp::Rcout << "[get_unik_particles]: tmp_unik_particles.size() = " << tmp_unik_particles.size() << endl;
  // while we're at it, let's sort the particle set
  arma::vec p_star_vec = arma::zeros<vec>(num_unik_particles);
  arma::uvec indices(num_unik_particles);
  for(int l = 0; l < num_unik_particles; l++){
    p_star_vec(l) = tmp_p_star[l];
  }
  indices = arma::sort_index(p_star_vec, "descend");
  for(int l = 0; l < num_unik_particles; l++){
    unik_particles.push_back(tmp_unik_particles[indices(l)]);
    p_star.push_back(tmp_p_star[indices(l)]);
    counts.push_back(tmp_counts[indices(l)]);
  }
  
}
*/

// function to compute entropy
// when we replace the (current_l)^th particle with the candidate_clusters
//
double Entropy(unsigned current_l, Partition* candidate_particle, std::vector<LPPartition> particle_set, std::vector<double> w){
  unsigned L = particle_set.size();
  // need to loop over to extract the unique partitions
  std::vector<LPPartition> unik_particles;
  std::vector<double> p_star;

  unik_particles.push_back(candidate_particle);
  p_star.push_back(w[current_l]);

  // in a sense, we are replacing particle_set[current_l] with candidate_particle
  // by adding it to the unik_particles vector

  int num_unik_particles = 1;
  int counter = 0;
  for(unsigned l = 0; l < L; l++){ // loop over current particle set
	counter = 0;
	if(l != current_l){
	//std::cout << "[getUnik]: l = " << l << std::endl;
	  for(unsigned ul = 0; ul < num_unik_particles; ul++){
	    if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
	    // l^th partition is equal to the ul^th unique partition
	    //std::cout << "particle " << l << " is equal to unik particle " << ul << std::endl;
	      p_star[ul] += w[l]; // update p_star
	      break;
        } else {
	      counter++;
        }
	  }

	  //std::cout << "[getUnik]: counter = " << counter << std::endl;
	  if(counter == num_unik_particles){
	    //std::cout << "we found a new unique particle!" << std::endl;
	    //particle_set[l]->Print_Partition();
	    // we have found a new unique particle
	    unik_particles.push_back(particle_set[l]);
	    p_star.push_back(w[l]);
	    num_unik_particles++;
      }
    }
  }
  double entropy = 0.0;
  //std::cout << "p_star = " ;
  for(unsigned ul = 0; ul < num_unik_particles; ul++){
	  //std::cout << p_star[ul] << " " ;
	  entropy += p_star[ul] * log(p_star[ul]);
  }
  //std::cout << std::endl;
  return -1.0 * entropy;

}

// a silly function to add up log-likelihoods and log-priors
//double total_log_post(LPPartition partition, const double a_sigma, const double nu_sigma){
double total_log_post(LPPartition partition, const double nu_sigma, const double lambda_sigma){
	//double log_post = total_log_like(partition, a_sigma, nu_sigma);
  double log_post = total_log_like(partition, nu_sigma, lambda_sigma);
  for(int k = 0; k < partition->K; k++){
    log_post += partition->log_prior[k];
  }
  
	//for(int k = 0; k < partition->K; k++){
	//	log_post += partition->log_like[k] + partition->log_prior[k];
	//}
	return log_post;
}

//double total_log_like(LPPartition partition, const double a_sigma, const double nu_sigma){
double total_log_like(LPPartition partition, const double nu_sigma, const double lambda_sigma){

  double log_like = 0.0;
  double log_det = 0.0;
  double quad_form = 0.0;
  
  for(int k = 0; k < partition->K; k++){
    log_det += partition->log_det_Omegay[k];
    quad_form += partition->y_Omegay_y[k];
  }
  
  //log_like = 0.5 * log_det - (a_sigma + partition->nObs/2) * log( (nu_sigma + quad_form)/2);
  log_like = 0.5 * log_det - ( (nu_sigma + partition->nObs)/2) * log( (nu_sigma * lambda_sigma + quad_form)/2);

  return log_like;
}
double total_log_prior(LPPartition partition){
  double log_prior = 0.0;
  for(int k = 0; k < partition->K; k++){
	log_prior += partition->log_prior[k];
  }
  return log_prior;
}

double alpha_bar_func(std::vector<int> new_cluster, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2)
{
  int n_k = new_cluster.size();
  if(n_k == 1){
    return(gamma_l->alpha_hat[new_cluster[0]]);
  } else{
    double tmp_alpha_bar = 0.0;
    for(int i = 0; i < n_k; i++){
      tmp_alpha_bar += (1 - rho) * gamma_l->alpha_hat[new_cluster[i]];
    }
    return( ((1.0/a1) * tmp_alpha_bar)/(1/a2 + n_k * (1-rho)/a1));
  }
}

//void update_w(std::vector<LPPartition> particle_set, std::vector<double> &w, const int L, const double a_sigma, const double nu_sigma, const double lambda)
void update_w(std::vector<LPPartition> particle_set, std::vector<double> &w, const int L, const double nu_sigma, const double lambda_sigma, const double lambda)
{
  
  double max_log_post = 0.0;
  double tmp_log_post = 0.0;
  double tmp_norm = 0.0;
  double tmp_p_star = 0.0;
  //Rcpp::Rcout << "[update_w]: Entering" << endl;

  // First we need to identify the unique particles
  std::vector<LPPartition> unik_particles;
  unik_particles.push_back(particle_set[0]);
  std::vector<int> particle_assignment(L,-1); // tells us to which unique particle each element of particle_ste corresponds
  particle_assignment[0] = 0;
  std::vector<int> particle_counts;
  
  std::vector<double> p_star;
  std::vector<double> log_post;
  //Rcpp::Rcout << "[update_w]: About to compute total_log_post(particle_set[0])" << endl;
  //Rcpp::Rcout << total_log_post(particle_set[0]) << endl;
  //particle_set[0]->Print_Partition();
  //log_post.push_back(total_log_post(particle_set[0], a_sigma, nu_sigma));
  log_post.push_back(total_log_post(particle_set[0], nu_sigma, lambda_sigma));
  p_star.push_back(0);
  //max_log_post = total_log_post(particle_set[0], a_sigma, nu_sigma);
  max_log_post = total_log_post(particle_set[0], nu_sigma, lambda_sigma);
  //Rcpp::Rcout << "[update_w]: Ready to find unik_particles" << endl;
  
  int num_unik_particles = 1;
  int counter = 0;

  for(int l = 1; l < L; l++){ // loop over the particle set
    counter = 0;
    for(int ul = 0; ul < num_unik_particles; ul++){
      if(Partition_Equal(particle_set[l], unik_particles[ul]) == 1){
        // l^th particle is equal to the ul^th unique partition
        particle_assignment[l] = ul;
        break;
      } else{
        counter++;
      }
    } // closes loop over the unique particles
    if(counter == num_unik_particles){
      // we have found a new unique particle
      unik_particles.push_back(particle_set[l]);
      particle_assignment[l] = num_unik_particles;
      p_star.push_back(0.0); // for now we populate p_star with 0's
      //tmp_log_post = total_log_post(particle_set[l], a_sigma, nu_sigma);
      tmp_log_post = total_log_post(particle_set[l], nu_sigma, lambda_sigma);
      log_post.push_back(tmp_log_post);
      if(tmp_log_post > max_log_post){
        max_log_post = tmp_log_post;
      }
      num_unik_particles++;
    }
  }
  //Rcpp::Rcout << "[update_w] : num_unik_particles = " << num_unik_particles << endl;

  particle_counts.clear();
  particle_counts.resize(num_unik_particles,0);
  for(int l = 0; l < L; l++){
    particle_counts[particle_assignment[l]]++;
  }
  for(int ul = 0; ul < num_unik_particles;ul++){
    tmp_log_post = log_post[ul] - max_log_post;
    tmp_p_star = exp(1/lambda * tmp_log_post);
    tmp_norm += tmp_p_star;
    p_star[ul] = tmp_p_star;
  }
  for(int ul = 0; ul < num_unik_particles;ul++){
    p_star[ul] /= tmp_norm;
  }
  for(int l = 0; l < L; l++){
    w[l] = p_star[particle_assignment[l]]/( (double) particle_counts[particle_assignment[l]]);
  }
}




// island move just removes the min and max from the cluster
void get_island(split_info &si,LPPartition gamma_l, const arma::mat &A_block)
{
  //Rcpp::Rcout << "[get_island]: Entering" << endl;
  si.num_splits = 0;
  si.split_k.clear();
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int n_k = 0;
  arma::vec alpha_hat_cluster = arma::zeros<vec>(1);
  arma::uvec indices(1);

  std::vector<int> k_star;
  
  std::vector<int> remain; // holds the indices of elements that remain in cluster split_k
  
  std::vector<std::vector<int> > new_clusters; // initially holds connected components of split_k after the split. then the last element is the island being created
  
  std::vector<std::vector<int> > tmp_new_clusters;
  arma::mat A_tmp = arma::zeros<mat>(1,1);
  int K = gamma_l->K;
  for(int k = 0; k < K; k++){
    n_k = gamma_l->cluster_config[k];
    //Rcpp::Rcout << "[get_island]:  k = " << k << " n_k = " << n_k << arma::endl;
    if(n_k > 1){
      alpha_hat_cluster.reset();
      alpha_hat_cluster.set_size(n_k);
      indices.reset();
      indices.set_size(n_k);
      for(int i = 0; i < n_k; i++){
        alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[k][i]];
      }
      indices = arma::sort_index(alpha_hat_cluster, "ascend"); //
      
      
      remain.clear();
      // remove the minimum element
      for(int i = 1; i < n_k; i++){
        remain.push_back(gamma_l->clusters[k][indices(i)]);
      }
      //Rcpp::Rcout << "min element = " << gamma_l->clusters[k][indices(0)];
      //Rcpp::Rcout << "max element = " << gamma_l->clusters[k][indices(n_k - 1)];
      //Rcpp::Rcout << "remain has size" << remain.size() << endl;
     // Rcpp::Rcout << "remain[0] = " << remain[0] << endl;
      //Rcpp::Rcout << "remain[n_k-1] = " << remain[n_k - 1] << endl;
      
      
      // find the connected components of remain
      tmp_new_clusters.clear();
      new_clusters.clear();
      new_clusters.push_back(std::vector<int>(1, remain[0])); // create the first component, containing first element of remain
      for(int i = 1; i < remain.size(); i++){
        tmp_new_clusters.clear();
        tmp_new_clusters.push_back(std::vector<int>(1, remain[i]));
        for(int new_k = 0; new_k < new_clusters.size(); new_k++){
          A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), new_clusters[new_k].size(), tmp_new_clusters[0], new_clusters[new_k]);
          if(any(vectorise(A_tmp) == 1.0)){
            for(int ii = 0; ii < new_clusters[new_k].size(); ii++){
              tmp_new_clusters[0].push_back(new_clusters[new_k][ii]);
            }
          } else{
            tmp_new_clusters.push_back(new_clusters[new_k]);
          }
          new_clusters[new_k].clear();
        } // closes loop over components of remain_clusters
        new_clusters.clear();
        for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
          new_clusters.push_back(tmp_new_clusters[new_k]);
          tmp_new_clusters[new_k].clear();
        }
        tmp_new_clusters.clear();
      }
      // now add the island to the back of new_clusters
      new_clusters.push_back(std::vector<int>(1, gamma_l->clusters[k][indices(0)]));
      k_star.clear();
      k_star.resize(new_clusters.size(), -1);
      
      
      si.split_k.push_back(k);
      si.new_clusters.push_back(new_clusters);
      si.nearest_neighbor.push_back(k_star);
      si.num_splits++;
      
      // now remove the max element
      remain.clear();
      for(int i = 0; i < n_k-1; i++){
        remain.push_back(gamma_l->clusters[k][indices(i)]);
      }
      // find the connected components of remain
      tmp_new_clusters.clear();
      new_clusters.clear();
      new_clusters.push_back(std::vector<int>(1, remain[0])); // create the first component, containing first element of remain
      for(int i = 1; i < remain.size(); i++){
        tmp_new_clusters.clear();
        tmp_new_clusters.push_back(std::vector<int>(1, remain[i]));
        for(int new_k = 0; new_k < new_clusters.size(); new_k++){
          A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), new_clusters[new_k].size(), tmp_new_clusters[0], new_clusters[new_k]);
          if(any(vectorise(A_tmp) == 1.0)){
            for(int ii = 0; ii < new_clusters[new_k].size(); ii++){
              tmp_new_clusters[0].push_back(new_clusters[new_k][ii]);
            }
          } else{
            tmp_new_clusters.push_back(new_clusters[new_k]);
          }
          new_clusters[new_k].clear();
        } // closes loop over components of remain_clusters
        new_clusters.clear();
        for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
          new_clusters.push_back(tmp_new_clusters[new_k]);
          tmp_new_clusters[new_k].clear();
        }
        tmp_new_clusters.clear();
      }
      new_clusters.push_back(std::vector<int>(1, gamma_l->clusters[k][indices(n_k-1)]));
      k_star.clear();
      k_star.resize(new_clusters.size(), -1);
      
      si.split_k.push_back(k);
      si.new_clusters.push_back(new_clusters);
      si.nearest_neighbor.push_back(k_star);
      si.num_splits++;

    }
  }
}


void get_local(split_info &si, LPPartition gamma_l, const arma::mat &A_block)
{
  //Rcpp::Rcout << "[get_local]: Entering" << endl;
  si.num_splits = 0;
  si.split_k.clear();
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int n_k = 0;
  arma::vec alpha_hat_cluster = arma::zeros<vec>(1);
  arma::uvec indices(1);
  
  std::vector<int> k_star;
  
  std::vector<int> remain; // holds the indices of elements that remain in cluster split_k
  
  std::vector<std::vector<int> > new_clusters; // initially holds connected components of split_k after the split. then the last element is the island being created
  
  std::vector<std::vector<int> > tmp_new_clusters;
  arma::mat A_tmp = arma::zeros<mat>(1,1);
  int K = gamma_l->K;

  std::vector<int> tmp_nn;
  std::vector<double> tmp_distance;
  arma::vec distance(1);
  arma::uvec dist_indices(1);
  
  
  for(int k = 0; k < K; k++){
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){
      alpha_hat_cluster.reset();
      alpha_hat_cluster.set_size(n_k);
      indices.reset();
      indices.set_size(n_k);
      for(int i = 0; i < n_k; i++){
        alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[k][i]];
      }
      indices = arma::sort_index(alpha_hat_cluster, "ascend");
      // loop over the clusters and try to remove each element of the cluster individually
      for(int i = 0; i < n_k; i++){
        remain.clear();
        // remove the minimum element
        for(int ii = 0; ii < n_k; ii++){
          // remember to skip over the element labelled i
          if(ii != i) remain.push_back(gamma_l->clusters[k][indices(ii)]);
        }
        // find the connected components of remain
        tmp_new_clusters.clear();
        new_clusters.clear();
        new_clusters.push_back(std::vector<int>(1, remain[0])); // create the first component, containing first element of remain
        for(int ii = 1; ii < remain.size(); ii++){
          tmp_new_clusters.clear();
          tmp_new_clusters.push_back(std::vector<int>(1, remain[ii]));
          for(int new_k = 0; new_k < new_clusters.size(); new_k++){
            A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), new_clusters[new_k].size(), tmp_new_clusters[0], new_clusters[new_k]);
            if(any(vectorise(A_tmp) == 1.0)){
              for(int ix = 0; ix < new_clusters[new_k].size(); ix++){
                tmp_new_clusters[0].push_back(new_clusters[new_k][ix]);
              }
            } else{
              tmp_new_clusters.push_back(new_clusters[new_k]);
            }
            new_clusters[new_k].clear();
          } // closes loop over components of remain_clusters
          new_clusters.clear();
          for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
            new_clusters.push_back(tmp_new_clusters[new_k]);
            tmp_new_clusters[new_k].clear();
          }
          tmp_new_clusters.clear();
        }
        // now add the island to the back of new_clusters
        k_star.clear();
        k_star.resize(new_clusters.size(), -1); // up to this point k_star has length only equal to the number of connected components remaining in cluster k
        
        new_clusters.push_back(std::vector<int>(1, gamma_l->clusters[k][indices(i)]));
        
        tmp_nn.clear();
        tmp_distance.clear();
        // need to figure out what the nearest neighbor of the element being re-allocated is
        for(int kk = 0; kk < gamma_l->K; kk++){ // loop over the existing clusters
          if(kk != k){
            A_tmp = Submatrix(A_block, 1, gamma_l->cluster_config[kk], std::vector<int>(gamma_l->clusters[k][indices(i)]), gamma_l->clusters[kk]);
            if(any(vectorise(A_tmp) == 1)){
              // cluster kk is a neighbor of the island being considered
              tmp_nn.push_back(kk);
            }
          }
        }
        if(tmp_nn.size() > 0){
          distance.reset();
          dist_indices.reset();
          distance.set_size(tmp_nn.size());
          dist_indices.set_size(tmp_nn.size());
          for(int ii = 0; ii < tmp_nn.size(); ii++){
            distance(ii) = tmp_distance[ii];
          }
          dist_indices = arma::sort_index(distance, "ascend");
          k_star.push_back(tmp_nn[dist_indices(0)]);
        } else{
          k_star.push_back(-1);
        }
        si.split_k.push_back(k);
        si.new_clusters.push_back(new_clusters);
        si.nearest_neighbor.push_back(k_star);
        si.num_splits++;
      } // closes loop over the elements of cluster k
    } // closes if statement checking if cluster size > 1
  } // closes for loop over the clusters
}

/*
void best_island(LPPartition candidate, const int current_l, split_info& si, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a, const double eta, const double lambda)
{
  
  if(si.num_splits == 0){
    // no valid island candidates.
    // set candidate = particle_set[current_l]
    candidate->Copy_Partition(particle_set[current_l]);
    //Rcpp::Rcout << "  [best_island]: no island candidates. Returning gamma_l" << endl;
  } else{
    LPPartition max_candidate = new Partition(particle_set[current_l]); // holds the running "best" island candidate
    double max_objective = w[current_l]*total_log_post(max_candidate) + lambda*Entropy(current_l, max_candidate, particle_set, w);
    LPPartition tmp_candidate = new Partition(particle_set[current_l]); // the running candidate
    double tmp_objective = 0.0;

    
    // si contains just about everything we need
    for(int ix = 0; ix < si.num_splits; ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[current_l]);
      tmp_candidate->Split_Merge(si.split_k[ix], si.new_clusters[ix], si.nearest_neighbor[ix], ybar, T, A_block, rho, a, eta);
      tmp_objective = w[current_l]*total_log_post(tmp_candidate) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
      
      //Rcpp::Rcout << "current island proposal is :" << endl;
      //tmp_candidate->Print_Partition();
      //Rcpp::Rcout << "tmp_objective = " << tmp_objective << "   max_objective = " << max_objective << endl;
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
      }
    }
    candidate->Copy_Partition(max_candidate);
    delete max_candidate;
    delete tmp_candidate;
    //delete[] new_cluster1;
    //delete[] new_cluster2;
    
  }
}
*/
void get_border(split_info &si, LPPartition gamma_l, const arma::mat &A_block)
{
  
  si.num_splits = 0;
  si.split_k.clear();
  si.nearest_neighbor.clear();
  si.new_clusters.clear();
  
  if(gamma_l->K > 1){
    arma::vec distance = arma::zeros<vec>(1); // holds distance of individual alpha estimates from the cluster mean
    arma::uvec indices(1);
    
    double alpha_bar = 0.0;
    arma::mat A_sub = A_block;
    std::vector<std::vector<int> > new_clusters; // holds connected components of cluster split during border move
    std::vector<std::vector<int> > tmp_new_clusters; // temporarily holds connected component of the cluster that is split during border move
    std::vector<int> k_star;
    std::vector<int> possible_border(1); // holds id of the potential border elements
    std::vector<int> possible_distance(1); // holds distance of potential border elements to cluster mean
    int border = 0; // the actual border
    int split_k = 0;
    int first_index = 0; // index of first element in cluster being split which is not the border
    arma::mat A_tmp;

    for(int k = 0; k < gamma_l->K; k++){
      
      possible_border.clear();
      possible_distance.clear();
      //for(int kk = k+1; kk < gamma_l->K; kk++){
      for(int kk = 0; kk < gamma_l->K; kk++){
        if( (gamma_l->cluster_config[kk] > 1) && (kk != k)){
          A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[kk], gamma_l->clusters[k], gamma_l->clusters[kk]);
          if(any(vectorise(A_sub) == 1.0)){ // cluster k and kk are adjacent
            for(int i = 0; i < gamma_l->cluster_config[kk]; i++){ // loop over element of cluster kk to see if it is a potential border
              if(any(A_sub.col(i) == 1.0)){
                possible_border.push_back(gamma_l->clusters[kk][i]);
                possible_distance.push_back(abs(alpha_bar - gamma_l->alpha_hat[gamma_l->clusters[kk][i]]));
              }
            }
          }
        }
      }
      if(possible_border.size() > 0){
        distance.reset();
        indices.reset();
        distance.set_size(possible_border.size());
        indices.set_size(possible_border.size());
        for(int i = 0; i < possible_border.size(); i++){
          distance(i) = possible_distance[i];
        }
        indices = arma::sort_index(distance, "ascend");
        border = possible_border[indices(0)]; // block group that would be merged with k.
        //Rcpp::Rcout << "[get_border]: blockgroup" << border << " will move from cluster " << gamma_l->cluster_assignment[border] << " to " << k << endl;
        split_k = gamma_l->cluster_assignment[border]; // this is the cluster we are removing border from
        //Rcpp::Rcout << "[get_border]: k = " << k << " split_k = " << split_k << " border = " << border << endl;
        for(int i = 0; i < gamma_l->cluster_config[split_k]; i++){
          if(gamma_l->clusters[split_k][i] != border){
            first_index = i;
            break;
          }
        }
        // gamma_l->clusters[split_k][first_index] is the first element of cluster being split that is not the border element
        new_clusters.clear();
        tmp_new_clusters.clear();
        new_clusters.push_back(std::vector<int>(1, gamma_l->clusters[split_k][first_index]));
        for(int i = 0; i < gamma_l->cluster_config[split_k]; i++){
          if( (gamma_l->clusters[split_k][i] != border) && (i != first_index)){
            tmp_new_clusters.clear();
            tmp_new_clusters.push_back(std::vector<int>(1, gamma_l->clusters[split_k][i]));
            for(int new_k = 0; new_k < new_clusters.size(); new_k++){
              A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), new_clusters[new_k].size(), tmp_new_clusters[0], new_clusters[new_k]);
              if(any(vectorise(A_tmp) == 1.0)){
                for(int ii = 0; ii < new_clusters[new_k].size(); ii++){
                  tmp_new_clusters[0].push_back(new_clusters[new_k][ii]);
                }
              } else{
                tmp_new_clusters.push_back(new_clusters[new_k]);
              }
              new_clusters[new_k].clear();
            }
            new_clusters.clear();
            for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
              new_clusters.push_back(tmp_new_clusters[new_k]);
              tmp_new_clusters[new_k].clear();
            }
          } // closes if statement that skips over the border
        }
        // new clusters is now done. It currently contains connected components of split_k after element border is removed
        k_star.clear();
        k_star.resize(new_clusters.size(),-1);
        // add border to end of new_cluster
        new_clusters.push_back(std::vector<int>(1, border));
        // add k to indicate that in fact border gets moved to cluster k
        k_star.push_back(k);
        if(new_clusters.size() != k_star.size()){
          Rcpp::Rcout << "[get_border]: k_star and new_cluster are of different sizes!" << endl;
        }
        
        
        //Rcpp::Rcout << "[get_border]: Proposing the following border move" << endl;
        //Rcpp::Rcout << "  cluster " << split_k << " split into : " << endl;
        //for(int new_k = 0; new_k < new_clusters.size(); new_k++){
        //  Rcpp::Rcout << "  new_cluster " << new_k << " of size " << new_clusters[new_k].size() << " with neighbor " << k_star[new_k] << " : " << endl;
        //  for(int i = 0; i < new_clusters[new_k].size(); i++){
        //    Rcpp::Rcout << new_clusters[new_k][i] << " " ;
        //  }
        //  Rcpp::Rcout << endl;
        //}
        si.new_clusters.push_back(new_clusters);
        si.nearest_neighbor.push_back(k_star);
        si.split_k.push_back(split_k);
        si.num_splits++;
      } // closes if statement checking if a border move is possible
    } // closes loop over the clusters
  } // closes if statement checking that gamma_l->K >1

}


/*
void get_merge(split_info& si, LPPartition gamma_l, const arma::mat &A_block)
{
  si.num_splits = 0;
  si.split_k.clear();
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  double alpha_bar = 0.0; // holds the cluster mean
  std::vector<int> tmp_neighbor(1); // holds potential neighbors that are being merged into k
  std::vector<double> tmp_distance(1); // holds distance from k to potential neighbors
  arma::vec distance(1);
  arma::uvec indices(1);
  arma::mat A_sub = arma::zeros<mat>(1,1);
  int split_k = 0;
  std::vector<std::vector<int> > new_clusters;
  
  if(gamma_l->K > 1){
    for(int k = 0; k < gamma_l->K; k++){
      alpha_bar = gamma_l->alpha_bar[k];
      tmp_neighbor.clear();
      tmp_distance.clear();
      for(int kk = k+1;kk < gamma_l->K; kk++){
        A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[kk], gamma_l->clusters[k], gamma_l->clusters[kk]);
        if(any(vectorise(A_sub) == 1)){ // clusters k and kk are neighbors
          tmp_neighbor.push_back(kk);
          tmp_distance.push_back(abs(alpha_bar - gamma_l->alpha_bar[kk]));
        }
      }
      if(tmp_neighbor.size() > 0){
        distance.reset();
        indices.reset();
        distance.set_size(tmp_neighbor.size());
        indices.set_size(tmp_neighbor.size());
        for(int i = 0; i < tmp_neighbor.size(); i++){
          distance(i) = tmp_distance[i];
        }
        indices = arma::sort_index(distance, "ascend");
        split_k = tmp_neighbor[indices(0)];
        
        new_clusters.clear();
        new_clusters.push_back(gamma_l->clusters[split_k]);
        
        si.split_k.push_back(split_k);
        si.new_clusters.push_back(new_clusters);
        si.nearest_neighbor.push_back(std::vector<int>(1,k));
        si.num_splits++;
      } // closes if statement checking that cluster k is adjacent to another cluster
    } // closes for loop over the clusters
  } // closes if statement checking that there are at least 2 clusters
}
*/

void get_merge(merge_info &mi, LPPartition gamma_l, const arma::mat &A_block)
{
  mi.num_merges = 0;
  mi.rec_k.clear();
  mi.donor_k.clear();
  
  double alpha_bar = 0.0; // holds the cluster mean
  std::vector<int> tmp_neighbor(1); // holds potential neighbors that are being merged into k
  std::vector<double> tmp_distance(1); // holds distance from k to potential neighbors
  arma::vec distance(1);
  arma::uvec indices(1);
  arma::mat A_sub = arma::zeros<mat>(1,1);
  std::vector<std::vector<int> > new_clusters;
  
  if(gamma_l->K > 1){
    for(int k = 0; k < gamma_l->K; k++){
      alpha_bar = gamma_l->alpha_bar[k];
      tmp_neighbor.clear();
      tmp_distance.clear();
      for(int kk = k+1;kk < gamma_l->K; kk++){
        A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[kk], gamma_l->clusters[k], gamma_l->clusters[kk]);
        if(any(vectorise(A_sub) == 1)){ // clusters k and kk are neighbors
          tmp_neighbor.push_back(kk);
          tmp_distance.push_back(abs(alpha_bar - gamma_l->alpha_bar[kk]));
        }
      }
      if(tmp_neighbor.size() > 0){
        distance.reset();
        indices.reset();
        distance.set_size(tmp_neighbor.size());
        indices.set_size(tmp_neighbor.size());
        for(int i = 0; i < tmp_neighbor.size(); i++){
          distance(i) = tmp_distance[i];
        }
        indices = arma::sort_index(distance, "ascend");
        
        mi.rec_k.push_back(k);
        mi.donor_k.push_back(tmp_neighbor[indices(0)]);
        mi.num_merges++;
      }
    }
  }
}

void get_spectral_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double split_frac){
  
  // re-set everything in si
  si.num_splits = 0;
  si.split_k.clear();
  for(int i = 0; i < si.new_clusters.size(); i++){
    si.new_clusters[i].clear();
    si.nearest_neighbor[i].clear();
  }
  
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int orig_K = gamma_l->K;
  int max_splits = 2;
  int split_k = 0;
  int n_k = 1;
  /*
   Initialize all of the stuff for spectral clustering here
   */
  arma::mat A_block_k = arma::zeros<mat>(n_k, n_k);
  arma::vec alpha_hat_cluster = arma::zeros<vec>(n_k);
  arma::mat alpha_sim = arma::zeros<mat>(n_k,n_k);
  arma::mat I_k = arma::zeros<mat>(n_k,n_k);
  I_k.eye();
  arma::mat alpha_dist = arma::zeros<mat>(n_k,n_k);
  arma::mat W_alpha_cl = arma::zeros<mat>(n_k,n_k);
  arma::mat Dinv_sqrt = arma::zeros<mat>(n_k,n_k);
  arma::mat L = arma::zeros<mat>(n_k,n_k);
  arma::mat eigvec = arma::zeros<mat>(n_k,n_k); // holds the eigenvectors of L
  arma::vec eigval = arma::zeros<vec>(n_k); // holds the eigenvalues of L
  arma::mat U = arma::zeros<mat>(n_k,2); // will hold the first several eigenvalues of L
  arma::mat means = arma::zeros<mat>(n_k,2); // the cluster means passed to kmeans
  bool status = true;
  arma::vec distance = arma::zeros<vec>(2);
  arma::uvec indices(2);

  std::vector<std::vector<int> > tmp_new_clusters; // holds the new clusters in unsorted order
  std::vector<std::vector<int> > new_clusters; // holds new clusters in sorted order.
  
  arma::vec tmp_nc_size = arma::zeros<vec>(1); // holds the size of the new clusters
  arma::uvec tmp_nc_size_indices(1); // used to sort tmp_nc_size
  arma::mat tmp_A = arma::zeros<mat>(1,1); // used to determine if new cluster is adjacent to an existing cluster
  double new_alpha_bar = 0.0;
  int nn_counter = 0; // counts number of neighbors of each new cluster
  std::vector<int> nn(1); // holds cluster labels of neighbors of each new cluster
  std::vector<double> nn_dist(1); // holds distance of each neighbor to new cluster
  arma::vec nn_dist_vec = arma::zeros<vec>(1); // arma vector version of nn_dist
  arma::uvec nn_dist_indices(1); // for soring nn_dist_vec
  std::vector<int> nearest_neighbor; // holds the cluster labels of the nearest neighbors
  
  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    //Rcpp::Rcout << "[get_spectral_splits]: Trying to split cluster " << k << endl;
    if(n_k > 1){
      if(2 >= ceil(split_frac * sqrt(n_k))) max_splits = 3;
      else max_splits = ceil(split_frac * sqrt(n_k));
      //Rcpp::Rcout << "max splits = " << max_splits << endl;
      for(int num_splits = 2; num_splits < max_splits; num_splits++){
        //Rcpp::Rcout << "    starting to split into " << num_splits << "  pieces" << endl;
        // re-size all of the things needed for the spectral clustering
        A_block_k.set_size(n_k,n_k);
        alpha_hat_cluster.set_size(n_k);
        alpha_sim.set_size(n_k,n_k);
        I_k.set_size(n_k,n_k);
        I_k.eye();
        alpha_dist.set_size(n_k,n_k);
        W_alpha_cl.set_size(n_k,n_k);
        Dinv_sqrt.set_size(n_k,n_k);
        L.set_size(n_k,n_k);
        eigvec.set_size(n_k,n_k);
        eigval.set_size(n_k);
        U.set_size(n_k,num_splits);
        means.set_size(n_k, num_splits); // I think this should be means.set_size(num_splits, num_splits)
        //means.set_size(num_splits, num_splits); // holds centroids found by k-means. we are using num_splits-dimensinoal vectors and finding num_splits many clusters
        status = true;
        distance.set_size(num_splits);
        indices.set_size(num_splits);
      
        tmp_new_clusters.clear(); // holds the new clusters is unsorted order
        tmp_new_clusters.resize(num_splits);
        
        tmp_nc_size.set_size(num_splits); // holds size of the new clusters
        tmp_nc_size_indices.set_size(num_splits); // for sorting the new clusters based on size
        
        new_clusters.clear(); // gets written to si. records new clusters sorted by size
        new_clusters.resize(num_splits);
        
        // stuff used to compute the nearest neighbor of each new cluster
        nn_dist.resize(num_splits);
        nn_dist_vec.set_size(num_splits);
        nn_dist_indices.set_size(num_splits);
        nearest_neighbor.resize(num_splits,-1); // holds nearest neighbors.

        /* Spectral Clustering Begins */
        A_block_k = Submatrix(A_block, n_k, n_k, gamma_l->clusters[split_k], gamma_l->clusters[split_k]);
        for(int i = 0; i < n_k; i++){
          alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[split_k][i]];
        }
        
        alpha_dist = Distance_matrix(alpha_hat_cluster,n_k);
        alpha_sim = exp(-1.0 * square(alpha_dist)/(2 * arma::var(alpha_hat_cluster)));
        W_alpha_cl = I_k + alpha_sim % A_block_k;
        Dinv_sqrt = arma::diagmat(1/sqrt(arma::sum(W_alpha_cl,1)));
        L = I_k - Dinv_sqrt * W_alpha_cl * Dinv_sqrt;
        arma::eig_sym(eigval, eigvec,L);
        U = eigvec.cols(0,num_splits-1);
        U = arma::diagmat(1/sqrt(arma::sum(arma::square(U),1))) * U;
        status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
        if(status == false){
          Rcpp::Rcout << "kmeans failed!!" << endl;
        } else{
          for(int i = 0; i < n_k; i++){
            distance.zeros();
            for(int k = 0; k < num_splits; k++){
              distance(k) = arma::norm(U.row(i).t() - means.col(k)); // measures distance from point i in the cluster to the k^th new sub-cluster
            }
            indices = sort_index(distance);
            tmp_new_clusters[indices(0)].push_back(gamma_l->clusters[split_k][i]);
          }
          // now that we have new_clusters
          // let us arrange them by size
          tmp_nc_size.set_size(num_splits);
          tmp_nc_size_indices.set_size(num_splits);
          for(int nc_ix = 0; nc_ix < num_splits; nc_ix++){
            tmp_nc_size(nc_ix) = tmp_new_clusters[nc_ix].size();
          }
          tmp_nc_size_indices = arma::sort_index(tmp_nc_size, "ascend");
          for(int nc_ix = 0; nc_ix < num_splits; nc_ix++){
            new_clusters[nc_ix] = tmp_new_clusters[tmp_nc_size_indices(nc_ix)];
            nearest_neighbor[nc_ix] = -1;
          }
          //Rcpp::Rcout << "    sorted new clusters by size!" << endl;
          /* Spectral Cluster Ends */
          /* now we can identify the nearest neighbor */
          if(orig_K > 1){
            for(int nc_ix = 0; nc_ix < new_clusters.size(); nc_ix++){ // loop over the new clusters
              nn_counter = 0;
              nn.clear(); // holds cluster labels of clusters adjacent to new cluster
              nn_dist.clear(); // holds distances from new clusters to its neighbors
              
              new_alpha_bar = alpha_bar_func(new_clusters[nc_ix], gamma_l, T, A_block, rho, a1, a2);
              for(int kk = 0; kk < orig_K; kk++){ // loop over the existing clusters
                if(kk != split_k){
                  tmp_A = Submatrix(A_block, new_clusters[nc_ix].size(), gamma_l->cluster_config[kk], new_clusters[nc_ix], gamma_l->clusters[kk]);
                  if(any(vectorise(tmp_A) == 1)){ // cluster kk is adjacent to new cluster nc_ix
                    nn.push_back(kk); // add kk to running list of neighbors of nc_ix
                    nn_dist.push_back(abs(new_alpha_bar - gamma_l->alpha_bar[kk])); // add distance to kk to running list
                    nn_counter++; // increment counter on the number of neighbors
                  }
                }
              }// closes loop over existing clusters
              if(nn_counter > 0){
                nn_dist_vec.set_size(nn_counter);
                nn_dist_indices.set_size(nn_counter);
                for(int i = 0; i < nn_counter; i++){
                  nn_dist_vec(i) = nn_dist[i];
                }
                nn_dist_indices = arma::sort_index(nn_dist_vec, "ascend"); // nn_indices(0) gives index, within nn and nn_dist* of nearest neighbor
                nearest_neighbor[nc_ix] = nn[nn_dist_indices(0)];
              } else{ // new_cluster nc_ix has no neighbors
                nearest_neighbor[nc_ix] = -1;
              }
            } // closes loop over the new clusters
            //Rcpp::Rcout << "    got the nearest neighbors!" << endl;
          } else{
            for(int nc_ix = 0; nc_ix < new_clusters.size(); nc_ix++){
              nearest_neighbor[nc_ix] = -1;
            }
          } // closes if/else checking whether orig_K = 1
          
          /* time to update si */
          si.num_splits++;
          si.split_k.push_back(split_k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(nearest_neighbor);
        } // closes if/else checking whether kmeans failed or not
      } // closes for loop over num_splits
    } // closes if statement checking whether cluster is a singleton
  } // closes loop over all clusters
}

void get_tail_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double split_frac)
{
  si.num_splits = 0;
  si.split_k.clear();
  for(int i = 0; i < si.new_clusters.size(); i++){
    si.new_clusters[i].clear();
    si.nearest_neighbor[i].clear();
  }
  
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
  
  int n_k = 0;
  arma::vec alpha_hat_cluster = arma::zeros<vec>(1);
  arma::uvec alpha_hat_indices(1);
  
  std::vector<int> left_tail;
  std::vector<int> right_tail;
  std::vector<int> center;
  //std::vector<int> comb_tail;
  
  std::vector<int> sub_left_tail; // for sub-vectors of left-tail
  std::vector<int> sub_right_tail; // for sub-vectors of the right-tail
  std::vector<int> remain; // contains all of the indices that are not split

  // to determine the connected components
  arma::mat A_tmp; // sub-matrix whose row corresponds to new element being removed and columns are sub-clusters
  
  std::vector<std::vector<int> > left_new_clusters; // holds the connected components of sub_left_tail
  std::vector<std::vector<int> > right_new_clusters; // holds connected components of sub_right_tail
  std::vector<std::vector<int> > remain_clusters;
  
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily used to build the connected components.
  
  double tmp_alpha_bar = 0.0; // temporarily hold mean of the alpha-hat's in components of sub_*_tail
  std::vector<double> tmp_dist; // holds distance from the new sub-clusters to adjacent existing clusters
  std::vector<double> tmp_nn; // holds indices of potential nearest neighbors for new sub-clusters
  arma::vec tmp_dist_vec = arma::zeros<vec>(1); // arma vector version of tmp_dist
  arma::uvec tmp_dist_indices(1); // for soring tmp_dist_vec
  
  std::vector<int> left_k_star;
  std::vector<int> right_k_star;
  
  for(int k = 0; k < gamma_l->K; k++){ // loop over the clusters
    //Rcpp::Rcout << "[get_new_split]: Starting to split cluster " << k << endl;
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){ // we can only split cluster if it has more than 1 element
      alpha_hat_cluster.clear();
      alpha_hat_indices.clear();
      
      alpha_hat_cluster.set_size(n_k); // holds the alpha_hat's within the cluster
      alpha_hat_indices.set_size(n_k);
      
      for(int i = 0; i < n_k; i++){
        alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[k][i]];
      }
      alpha_hat_indices = arma::sort_index(alpha_hat_cluster, "ascend"); // sort in ascending order
      
      left_tail.clear();
      right_tail.clear();
      center.clear(); // middle portion of the distribution
      
      for(int i = 0; i < ceil(n_k*split_frac); i++){
        left_tail.push_back(gamma_l->clusters[k][alpha_hat_indices(i)]);
        right_tail.push_back(gamma_l->clusters[k][alpha_hat_indices(n_k-1 - i)]);
      }
      // figure out what is in the center of the distribution
      for(int i = ceil(n_k*split_frac); i < n_k-ceil(n_k*split_frac); i++){
        center.push_back(gamma_l->clusters[k][alpha_hat_indices(i)]);
      }
      /*
      Rcpp::Rcout << "left-tail has size " << left_tail.size() << " : " ;
      for(int i = 0; i < left_tail.size(); i++) Rcpp::Rcout << left_tail[i] << " ";
      Rcpp::Rcout << endl;
      Rcpp::Rcout << "center has size " << center.size() << " : " ;
      for(int i = 0; i < center.size(); i++) Rcpp::Rcout << center[i] << " ";
      Rcpp::Rcout << endl;
      Rcpp::Rcout << "right tail has size " << right_tail.size() << " : " ;
      for(int i = 0; i < right_tail.size(); i++) Rcpp::Rcout << right_tail[i] << " " ;
      Rcpp::Rcout << endl;
      */
      
      // now try to remove just the left-tail
      sub_left_tail.clear();
      sub_left_tail.push_back(left_tail[0]);
      left_new_clusters.clear();
      left_new_clusters.push_back(std::vector<int>(1, left_tail[0])); // start left_new_clusters off
      for(int i = 1; i < left_tail.size(); i++){
        sub_left_tail.push_back(left_tail[i]);
        remain.clear();
        // first the remaining things from the left-tail
        for(int ii = i+1; ii < left_tail.size(); ii++){
          remain.push_back(left_tail[ii]);
        }
        // now the things from the center
        for(int ii = 0; ii < center.size(); ii++){
          remain.push_back(center[ii]);
        }
        // also everything in the right tail will remain
        for(int ii = 0; ii < right_tail.size(); ii++){
          remain.push_back(right_tail[ii]);
        }
        /*
         // print statements to check our progress
        Rcpp::Rcout << "  sub_left_tail has size " << sub_left_tail.size() << " : " ;
        for(int ii = 0; ii < sub_left_tail.size(); ii++){
          Rcpp::Rcout << sub_left_tail[ii] << " " ;
        }
        Rcpp::Rcout << endl;
        Rcpp::Rcout << "  remain has size " << remain.size() << " : ";
        for(int ii = 0; ii < remain.size(); ii++){
          Rcpp::Rcout << remain[ii] << " ";
        }
        Rcpp::Rcout << endl;
        */
        // now find the connected components of remain
        tmp_new_clusters.clear();
        remain_clusters.clear();
        remain_clusters.push_back(std::vector<int>(1,remain[0]));
        //Rcpp::Rcout << "    attempting to find connected components of remain" << endl;
        for(int ii = 1; ii < remain.size(); ii++){
          // start building the components again, this time with remain[ii]
          tmp_new_clusters.clear();
          tmp_new_clusters.push_back(std::vector<int>(1, remain[ii]));
          // loop over the existing connected components
          for(int new_k = 0; new_k < remain_clusters.size(); new_k++){
            //Rcpp::Rcout << "    new_k = " << new_k ;
            //A_tmp = Submatrix(A_block, 1, remain_clusters[new_k].size(), std::vector<int>(1,remain[ii]), remain_clusters[new_k]); // this should work but we'll use following to be safe
            A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), remain_clusters[new_k].size(), tmp_new_clusters[0], remain_clusters[new_k]);
            if(any(vectorise(A_tmp) == 1.0)){
              // something in remain_clusters[new_k] is adjacent to remain[ii] so we should combine these clusters
              //Rcpp::Rcout << "   connected to tmp_new_clusters[0]!" << endl;
              for(int ix = 0; ix < remain_clusters[new_k].size(); ix++){
                tmp_new_clusters[0].push_back(remain_clusters[new_k][ix]);
              }
            } else{ // remain_clusters[new_k] remains its own distinct component
              //Rcpp::Rcout << "    not connected to tmp_new_clusters[0]!" << endl;
              tmp_new_clusters.push_back(remain_clusters[new_k]);
            }
            remain_clusters[new_k].clear();
          } // closes loop over elements of remain_clusters
          // update remain_clusters
          remain_clusters.clear();
          for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
            remain_clusters.push_back(tmp_new_clusters[new_k]);
            tmp_new_clusters[new_k].clear();
          }
          tmp_new_clusters.clear();
        } // closes loop over elements of remain used to determine connected components of remain
        /*
        // print statements to check our progress
        Rcpp::Rcout << "Connected components of remain : " << endl;
        for(int new_k = 0; new_k < remain_clusters.size(); new_k++){
          Rcpp::Rcout << "   component " << new_k << " with size " << remain_clusters[new_k].size() << " : ";
          for(int ii = 0; ii < remain_clusters[new_k].size(); ii++){
            Rcpp::Rcout << remain_clusters[new_k][ii] << " ";
          }
          Rcpp::Rcout << endl;
        }
        */
        
        // now it is time to find connected components of sub_left_tail
        tmp_new_clusters.clear();
        tmp_new_clusters.push_back(std::vector<int>(1,left_tail[i]));
        for(int new_k = 0; new_k < left_new_clusters.size(); new_k++){
          /*
          // print statements to check our progress
          //Rcpp::Rcout << "left_new_clusters[new_k] is: " ;
          //for(int ii = 0; ii < left_new_clusters[new_k].size(); ii++) Rcpp::Rcout << left_new_clusters[new_k][ii] << " ";
          //Rcpp::Rcout << endl;
          //Rcpp::Rcout << "tmp_new_clusters[0] is : " ;
          //for(int ii = 0; ii < tmp_new_clusters[0].size(); ii++) Rcpp::Rcout << tmp_new_clusters[0][ii] << " ";
          //Rcpp::Rcout << endl;
          */
          A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), left_new_clusters[new_k].size(), tmp_new_clusters[0], left_new_clusters[new_k]);
          if(any(vectorise(A_tmp) == 1.0)){
            for(int ii = 0; ii < left_new_clusters[new_k].size(); ii++){
              tmp_new_clusters[0].push_back(left_new_clusters[new_k][ii]);
            }
          } else{
            tmp_new_clusters.push_back(left_new_clusters[new_k]);
          }
          left_new_clusters[new_k].clear();
        }
        
        left_new_clusters.clear();
        left_k_star.clear();
        for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
          left_new_clusters.push_back(tmp_new_clusters[new_k]);
          // now figure out the connected components
          tmp_alpha_bar = alpha_bar_func(left_new_clusters[new_k], gamma_l, T, A_block, rho, a1, a2);
          tmp_nn.clear();
          tmp_dist.clear();
          for(int kk = 0; kk < gamma_l->K; kk++){
            if(kk != k){
              A_tmp = Submatrix(A_block, left_new_clusters[new_k].size(), gamma_l->cluster_config[kk], left_new_clusters[new_k], gamma_l->clusters[kk]);
              if(any(vectorise(A_tmp == 1.0))){
                tmp_nn.push_back(kk);
                tmp_dist.push_back(abs(tmp_alpha_bar - gamma_l->alpha_bar[kk]));
              }
            }
          }
          if(tmp_nn.size() > 0){
            tmp_dist_vec.reset();
            tmp_dist_indices.reset();
            tmp_dist_vec.set_size(tmp_dist.size());
            tmp_dist_indices.set_size(tmp_dist.size());
            for(int kk = 0; kk < tmp_dist.size(); kk++){
              tmp_dist_vec(kk) = tmp_dist[kk];
            }
            tmp_dist_indices = sort_index(tmp_dist_vec,"ascend");
            left_k_star.push_back(tmp_nn[tmp_dist_indices(0)]);
          } else{
            left_k_star.push_back(-1);
          }
        } // closes loop that adds components of sub_left_tail to left_new_clusters and finds their nearest neighbors
        /*
        // print statements used to check progress
        Rcpp::Rcout << "left_new_clusters:" << endl;
        for(int new_k = 0; new_k < left_new_clusters.size(); new_k++){
          Rcpp::Rcout << "    new_cluster " << new_k << " with size " << left_new_clusters[new_k].size() << " and neighbor " << left_k_star[new_k] << ":" << endl;
          for(int ii = 0; ii < left_new_clusters[new_k].size(); ii++){
            Rcpp::Rcout << left_new_clusters[new_k][ii] << " ";
          }
          Rcpp::Rcout << endl;
        }
        */
        // now we need to add things to si
        si.split_k.push_back(k);
        si.new_clusters.push_back(left_new_clusters);
        si.nearest_neighbor.push_back(left_k_star);
        // now add the elements from remaining portion
        for(int new_k = 0; new_k < remain_clusters.size(); new_k++){
          si.new_clusters[si.new_clusters.size()-1].push_back(remain_clusters[new_k]);
          si.nearest_neighbor[si.nearest_neighbor.size()-1].push_back(-1);
        }
        si.num_splits++;
      } // closes loop over elements of left_tail
      // now try to remove just the right tail
      sub_right_tail.clear();
      sub_right_tail.push_back(right_tail[0]);
      right_new_clusters.clear();
      right_new_clusters.push_back(std::vector<int>(1, right_tail[0])); // start right_new_clusters off
      for(int i = 1; i < right_tail.size(); i++){
        sub_right_tail.push_back(right_tail[i]);
        remain.clear();
        // first the remaining things from the right-tail
        for(int ii = i+1; ii < right_tail.size(); ii++){
          remain.push_back(right_tail[ii]);
        }
        // now the things from the center
        for(int ii = 0; ii < center.size(); ii++){
          remain.push_back(center[ii]);
        }
        // also everything in the left tail will remain
        for(int ii = 0; ii < left_tail.size(); ii++){
          remain.push_back(left_tail[ii]);
        }
        /*
        // print statement to check progress
        Rcpp::Rcout << "  sub_right_tail has size " << sub_right_tail.size() << " : " ;
        for(int ii = 0; ii < sub_right_tail.size(); ii++){
          Rcpp::Rcout << sub_right_tail[ii] << " " ;
        }
        Rcpp::Rcout << endl;
        Rcpp::Rcout << "  remain has size " << remain.size() << " : ";
        for(int ii = 0; ii < remain.size(); ii++){
          Rcpp::Rcout << remain[ii] << " ";
        }
        Rcpp::Rcout << endl;
        */
        // now find the connected components of remain
        tmp_new_clusters.clear();
        remain_clusters.clear();
        remain_clusters.push_back(std::vector<int>(1,remain[0]));
        //Rcpp::Rcout << "    attempting to find connected components of remain" << endl;
        for(int ii = 1; ii < remain.size(); ii++){
          // start building the components again, this time with remain[ii]
          tmp_new_clusters.clear();
          tmp_new_clusters.push_back(std::vector<int>(1, remain[ii]));
          // loop over the existing connected components
          for(int new_k = 0; new_k < remain_clusters.size(); new_k++){
            //Rcpp::Rcout << "    new_k = " << new_k ;
            //A_tmp = Submatrix(A_block, 1, remain_clusters[new_k].size(), std::vector<int>(1,remain[ii]), remain_clusters[new_k]); // this should work but we'll use following to be safe
            A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), remain_clusters[new_k].size(), tmp_new_clusters[0], remain_clusters[new_k]);
            if(any(vectorise(A_tmp) == 1.0)){
              // something in remain_clusters[new_k] is adjacent to remain[ii] so we should combine these clusters
              //Rcpp::Rcout << "   connected to tmp_new_clusters[0]!" << endl;
              for(int ix = 0; ix < remain_clusters[new_k].size(); ix++){
                tmp_new_clusters[0].push_back(remain_clusters[new_k][ix]);
              }
            } else{ // remain_clusters[new_k] remains its own distinct component
                    //Rcpp::Rcout << "    not connected to tmp_new_clusters[0]!" << endl;
              tmp_new_clusters.push_back(remain_clusters[new_k]);
            }
            remain_clusters[new_k].clear();
          } // closes loop over elements of remain_clusters
            // update remain_clusters
          remain_clusters.clear();
          for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
            remain_clusters.push_back(tmp_new_clusters[new_k]);
            tmp_new_clusters[new_k].clear();
          }
          tmp_new_clusters.clear();
        } // closes loop over elements of remain used to determine connected components of remain
        /*
        // print statements for checking progress
        //Rcpp::Rcout << "Connected components of remain : " << endl;
        //for(int new_k = 0; new_k < remain_clusters.size(); new_k++){
        //  Rcpp::Rcout << "   component " << new_k << " with size " << remain_clusters[new_k].size() << " : ";
        //   for(int ii = 0; ii < remain_clusters[new_k].size(); ii++){
        //    Rcpp::Rcout << remain_clusters[new_k][ii] << " ";
        //   }
        //   Rcpp::Rcout << endl;
        // }
        */
        // now it is time to find connected components of sub_right_tail
        tmp_new_clusters.clear();
        tmp_new_clusters.push_back(std::vector<int>(1,right_tail[i]));
        for(int new_k = 0; new_k < right_new_clusters.size(); new_k++){
          /*
          // print statements for checking progress
          Rcpp::Rcout << "right_new_clusters[new_k] is: " ;
          for(int ii = 0; ii < right_new_clusters[new_k].size(); ii++) Rcpp::Rcout << right_new_clusters[new_k][ii] << " ";
          Rcpp::Rcout << endl;
          Rcpp::Rcout << "tmp_new_clusters[0] is : " ;
          for(int ii = 0; ii < tmp_new_clusters[0].size(); ii++) Rcpp::Rcout << tmp_new_clusters[0][ii] << " ";
          Rcpp::Rcout << endl;
          */
          A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), right_new_clusters[new_k].size(), tmp_new_clusters[0], right_new_clusters[new_k]);
          if(any(vectorise(A_tmp) == 1.0)){
            for(int ii = 0; ii < right_new_clusters[new_k].size(); ii++){
              tmp_new_clusters[0].push_back(right_new_clusters[new_k][ii]);
            }
          } else{
            tmp_new_clusters.push_back(right_new_clusters[new_k]);
          }
          right_new_clusters[new_k].clear();
        }
        
        right_new_clusters.clear();
        right_k_star.clear();
        for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
          right_new_clusters.push_back(tmp_new_clusters[new_k]);
          // now figure out the connected components
          tmp_alpha_bar = alpha_bar_func(right_new_clusters[new_k], gamma_l, T, A_block, rho, a1, a2);
          tmp_nn.clear();
          tmp_dist.clear();
          for(int kk = 0; kk < gamma_l->K; kk++){
            if(kk != k){
              A_tmp = Submatrix(A_block, right_new_clusters[new_k].size(), gamma_l->cluster_config[kk], right_new_clusters[new_k], gamma_l->clusters[kk]);
              if(any(vectorise(A_tmp == 1.0))){
                tmp_nn.push_back(kk);
                tmp_dist.push_back(abs(tmp_alpha_bar - gamma_l->alpha_bar[kk]));
              }
            }
          }
          if(tmp_nn.size() > 0){
            tmp_dist_vec.reset();
            tmp_dist_indices.reset();
            tmp_dist_vec.set_size(tmp_dist.size());
            tmp_dist_indices.set_size(tmp_dist.size());
            for(int kk = 0; kk < tmp_dist.size(); kk++){
              tmp_dist_vec(kk) = tmp_dist[kk];
            }
            tmp_dist_indices = sort_index(tmp_dist_vec,"ascend");
            right_k_star.push_back(tmp_nn[tmp_dist_indices(0)]);
          } else{
            right_k_star.push_back(-1);
          }
        } // closes loop that adds components of sub_right_tail to right_new_clusters and finds their nearest neighbors
          /*
          // print statements for tracking progress
          Rcpp::Rcout << "right_new_clusters:" << endl;
          for(int new_k = 0; new_k < right_new_clusters.size(); new_k++){
            Rcpp::Rcout << "    new_cluster " << new_k << " with size " << right_new_clusters[new_k].size() << " and neighbor " << right_k_star[new_k] << ":" << endl;
            for(int ii = 0; ii < right_new_clusters[new_k].size(); ii++){
              Rcpp::Rcout << right_new_clusters[new_k][ii] << " ";
            }
            Rcpp::Rcout << endl;
          }
          */
        // now we need to add things to si
        si.split_k.push_back(k);
        si.new_clusters.push_back(right_new_clusters);
        si.nearest_neighbor.push_back(right_k_star);
        // now add the elements from remaining portion
        for(int new_k = 0; new_k < remain_clusters.size(); new_k++){
          si.new_clusters[si.new_clusters.size()-1].push_back(remain_clusters[new_k]);
          si.nearest_neighbor[si.nearest_neighbor.size()-1].push_back(-1);
        }
        si.num_splits++;
      }
      
    } // closes if checking whether n_k > 1
  } // closes loop over the clusters

  
}

void get_km_split(split_info &si, LPPartition gamma_l, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double split_frac)
{
    // re-set everything in si
  si.num_splits = 0;
  si.split_k.clear();
  for(int i = 0; i < si.new_clusters.size(); i++){
    si.new_clusters[i].clear();
    si.nearest_neighbor[i].clear();
  }
    
  si.new_clusters.clear();
  si.nearest_neighbor.clear();
    
  int orig_K = gamma_l->K;
  int max_splits = 2;
  int split_k = 0;
  int n_k = 1;
    
  // initialize the stuff needed for k-means clustering
  arma::mat U = arma::zeros<mat>(n_k, 1); // holds the data passed to k-means
  arma::mat means = arma::zeros<mat>(n_k,1); // holds data passed to k-means
  bool status = true; // tells us whether k-means was successful
    
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
    
  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
  std::vector<std::vector<int> > tmp_connected_components; // used to temporarily hold the connected components
    
    
  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
    
  arma::vec tmp_nc_size = arma::zeros<vec>(1); // holds size of the new clusters
  arma::uvec tmp_nc_size_indices(1); // used to sort tmp_nc_size
  arma::mat A_tmp = arma::zeros<mat>(1,1); // used to find nearest neighbors of newly formed clusters
  double tmp_alpha_bar = 0.0; // cluster mean for newly formed cluster
  std::vector<int> tmp_nn; // holds potential nearest neighbors
  std::vector<double> tmp_dist; // holds distance to potential nearest neighbors
  arma::vec tmp_dist_vec = arma::zeros<vec>(1);
  arma::uvec tmp_dist_indices(1);
  
  
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
    

  for(int k = 0; k < orig_K; k++){
    //Rcpp::Rcout << "[get_km_split]: k = " << k << endl;
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    if(n_k > 1){
      if(2 >= ceil(split_frac * sqrt(n_k))) max_splits = 3;
      else max_splits = ceil(sqrt(n_k) * split_frac);
      //Rcpp::Rcout << "max_splits = " << max_splits << endl;
      for(int num_splits = 2; num_splits <= max_splits; num_splits++){
        U.set_size(n_k, 1);
        means.set_size(1, num_splits);
        status = true;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        //Rcpp::Rcout << "    set U, means, status " << endl;
        
        // do k-means
        for(int ii = 0; ii < n_k; ii++){
          U(ii,0) = gamma_l->alpha_hat[gamma_l->clusters[k][ii]];
        }
        //Rcpp::Rcout << "Ready to start k-means" << endl;
        
        status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
  
        if(status == false){
          Rcpp::Rcout << "kmeans failed!!" << endl;
        } else{
          init_new_clusters.clear();
          init_new_clusters.resize(num_splits);
          //Rcpp::Rcout << "  init_new_clusters size = " << init_new_clusters.size() << endl;
          
          for(int i = 0; i < n_k; i++){
            cluster_dist.zeros();
            for(int new_k = 0; new_k < num_splits; new_k++){
              cluster_dist(new_k) = arma::norm(U.row(i).t() - means.col(new_k)); // distance from newly created cluster kk to point i
            }
            cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
            //Rcpp::Rcout << "i = " << i << "  goes to new cluster " << cluster_dist_indices(0) << endl;
            init_new_clusters[cluster_dist_indices(0)].push_back(gamma_l->clusters[split_k][i]);
          }
          

          // now loop over init_new_clusters and find the connected components
          tmp_new_clusters.clear();
          for(int kk = 0; kk < init_new_clusters.size(); kk++){
            connected_components.clear();
            connected_components.push_back(std::vector<int>(1, init_new_clusters[kk][0])); // start connected_components with first element of the newly found cluster
            tmp_connected_components.clear();
            
            for(int ii = 1; ii < init_new_clusters[kk].size(); ii++){
              tmp_connected_components.clear();
              tmp_connected_components.push_back(std::vector<int>(1, init_new_clusters[kk][ii]));
              for(int new_k = 0; new_k < connected_components.size(); new_k++){
                A_tmp = Submatrix(A_block, tmp_connected_components[0].size(), connected_components[new_k].size(), tmp_connected_components[0], connected_components[new_k]);
                if(any(vectorise(A_tmp) == 1)){
                  // something in connected_components[new_k] is adjacent to tmp_connected_components[0]
                  for(int ix = 0; ix < connected_components[new_k].size(); ix++){
                    tmp_connected_components[0].push_back(connected_components[new_k][ix]);
                  }
                } else{
                  tmp_connected_components.push_back(connected_components[new_k]);
                }
                connected_components[new_k].clear();
              } // closes loop over elements of connected components
              connected_components.clear();
              for(int new_k = 0; new_k < tmp_connected_components.size(); new_k++){
                connected_components.push_back(tmp_connected_components[new_k]);
                tmp_connected_components[new_k].clear();
              }
              tmp_connected_components.clear();
            } // closes loop over the elements of init_new_clusters[kk]
            // now that we have connected components of the sub-cluster discovered by k-means we can push it back to tmp_new_clusters
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the number of clusters found by k-means
          //Rcpp::Rcout << "   got connected components following k-means split " << endl;
          
          // at this point we have tmp_new_clusters. we'd like to sort them by size
          tmp_nc_size.reset();
          tmp_nc_size.set_size(tmp_new_clusters.size());
          tmp_nc_size_indices.reset();
          tmp_nc_size_indices.set_size(tmp_new_clusters.size());
          for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
            tmp_nc_size(new_k) = tmp_new_clusters[new_k].size();
          }
          tmp_nc_size_indices = arma::sort_index(tmp_nc_size, "ascend");
          new_clusters.clear();
          for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
            new_clusters.push_back(tmp_new_clusters[tmp_nc_size_indices(new_k)]);
          }
          //Rcpp::Rcout << "  number of new clusters = " << new_clusters.size() << endl;
          
          
          // now we can find the nearest neighbors
          k_star.clear();
          k_star.resize(new_clusters.size(), -1);
          for(int new_k = 0; new_k < new_clusters.size(); new_k++){
            tmp_alpha_bar = alpha_bar_func(new_clusters[new_k], gamma_l, T, A_block, rho, a1, a2);
            tmp_nn.clear();
            tmp_dist.clear();
            for(int kk = 0; kk < gamma_l->K; kk++){
              if(kk!=k){
                A_tmp = Submatrix(A_block, new_clusters[new_k].size(), gamma_l->cluster_config[kk], new_clusters[new_k], gamma_l->clusters[kk]);
                if(any(vectorise(A_tmp) == 1)){
                  tmp_nn.push_back(kk);
                  tmp_dist.push_back(abs(tmp_alpha_bar - gamma_l->alpha_bar[kk]));
                }
              }
            }
            if(tmp_nn.size() > 0){
              tmp_dist_vec.reset();
              tmp_dist_indices.reset();
              tmp_dist_vec.set_size(tmp_dist.size());
              tmp_dist_indices.set_size(tmp_dist.size());
              for(int kk = 0; kk < tmp_dist.size(); kk++){
                tmp_dist_vec(kk) = tmp_dist[kk];
              }
              tmp_dist_indices = arma::sort_index(tmp_dist_vec, "ascend");
              k_star[new_k] = tmp_nn[tmp_dist_indices(0)];
            } else{
              k_star[new_k] = -1;
            }
          }
          
          si.split_k.push_back(split_k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
        } // closes if/else checking that k-means succeeded
      } // closes loop over number of potential splits
    } // closes if statement checking that cluster can be split
    
    
    
  } // closes loop over the clusters of gamma_l

}



//void best_split(split_info &si, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double a_sigma, const double nu_sigma, const double eta, const double lambda)
void best_split(split_info &si, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda)
{
  
  if(si.num_splits > 0){
    //Rcpp::Rcout << "[best_split]: Will attempt " << si.num_splits << " different splits" << endl;
    LPPartition max_candidate = new Partition(particle_set[current_l]);
    //double max_objective = w[current_l]*total_log_post(max_candidate, a_sigma, nu_sigma) + lambda*Entropy(current_l, max_candidate, particle_set, w);
    double max_objective = w[current_l]*total_log_post(max_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, max_candidate, particle_set, w);
    LPPartition tmp_candidate = new Partition(particle_set[current_l]); // the running candidate
    double tmp_objective = 0.0;
    
    int break_flag = 0;
    
    std::vector<int> k_star; // holds the indices of nearest neighbors
    int num_new_clusters = 0;
    int split_k = 0;
    
    bool sanity_flag = true;
    
    for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
      //Rcpp::Rcout << "    [best_split]: split_ix = " << split_ix << endl;

      split_k = si.split_k[split_ix]; // which cluster is being split
      num_new_clusters = si.new_clusters[split_ix].size();
      
      k_star.clear();
      k_star.resize(num_new_clusters, -1); // note: simply re-sizing will not overwrite the values

 // print out some of the progress to make sure the split is valid
      /*
      Rcpp::Rcout << "Attempting to split cluster " << split_k << "into:" << endl;
      for(int i = 0; i < num_new_clusters;i++){
        Rcpp::Rcout << "  new cluster " << i << " of size " << si.new_clusters[split_ix][i].size() << " with neighbor " << k_star[i] << endl;
        for(int ii = 0; ii < si.new_clusters[split_ix][i].size(); ii++){
          Rcpp::Rcout << si.new_clusters[split_ix][i][ii] << " ";
        }
        Rcpp::Rcout << endl;
      }
*/
      // delete tmp_candidate and re-create it
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[current_l]);
      tmp_candidate->Split_Merge(split_k, si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
      sanity_flag = sanity_check(tmp_candidate);
      if(sanity_flag == false){
        Rcpp::Rcout << "[best_split]: something is wrong about this split!" << endl;
        Rcpp::Rcout << "   Attempted to split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " parts " << endl;
        for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
          Rcpp::Rcout << "    new cluster " << nc_ix << " of size " << si.new_clusters[split_ix][nc_ix].size() << " and neighbor " << si.nearest_neighbor[split_ix][nc_ix] << " : " << endl;
          for(int ii = 0; ii < si.new_clusters[split_ix][nc_ix].size(); ii++){
            Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][ii] << " " ;
          }
          Rcpp::Rcout << endl;
        }
        Rcpp::Rcout << "Before split candidate is" << endl;
        //particle_set[current_l]->Print_Partition(a_sigma, nu_sigma);
        particle_set[current_l]->Print_Partition(nu_sigma, lambda_sigma);
        Rcpp::Rcout << "The candidate is now" << endl;
        //tmp_candidate->Print_Partition(a_sigma, nu_sigma);
        tmp_candidate->Print_Partition(nu_sigma, lambda_sigma);
        Rcpp::stop("Terminating");
      }
      //tmp_objective = w[current_l]*total_log_post(tmp_candidate, a_sigma, nu_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
      tmp_objective = w[current_l]*total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);

      //Rcpp::Rcout << "tmp_objective = " << tmp_objective << "  max_objective = " << max_objective << endl;
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        //Rcpp::Rcout << "  max candidate updated. max_obj = " << max_objective  << endl;
        break_flag = 1;
      }
      // now we consider a sequence of progressive merges
/*
      if(break_flag == 0){
        for(int nc_ix = 0; nc_ix < num_new_clusters; nc_ix++){
          if(si.nearest_neighbor[split_ix][nc_ix] != -1){
            k_star[nc_ix] = si.nearest_neighbor[split_ix][nc_ix]; // as the loop progresses, we change one element of k_star at a time
            delete tmp_candidate;
            tmp_candidate = new Partition(particle_set[current_l]);
            tmp_candidate->Split_Merge(split_k, si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
            tmp_objective = w[current_l]*total_log_post(tmp_candidate, a_sigma, nu_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
            //Rcpp::Rcout << "    tmp_objective = " << tmp_objective << "  max_objective = " << max_objective << endl;
            if(tmp_objective > max_objective){
              delete max_candidate;
              max_candidate = new Partition(tmp_candidate);
              max_objective = tmp_objective;
              //Rcpp::Rcout << "  max candidate updated. max_obj = " << max_objective  << endl;

              break_flag = 1;
              break;
            }
          } // closes if statement checking if nearest neighbor of the next element is -1 or not; if it is -1 we don't try the progressive merge
          // not sure this is strictly necessary but we'll do it anyway
          if(break_flag == 1) break;
        }
      } // closes if statement that checks whether we should go into the progressive merges
*/
      if(break_flag == 1) break; // this will break us out of the loop over different proposed splits
    } // closes loop over the different proposed splits
    //Rcpp::Rcout << "Is max_candidate == particle_set[current_l]? " << Partition_Equal(particle_set[current_l], max_candidate) << endl;
    candidate->Copy_Partition(max_candidate);
    delete tmp_candidate;
    delete max_candidate;
  } else{
    candidate->Copy_Partition(particle_set[current_l]);
  } // closes if/else checking that there are splits proposed
}
//void best_merge(merge_info &mi, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double a_sigma, const double nu_sigma, const double eta, const double lambda)
void best_merge(merge_info &mi, LPPartition candidate, const int current_l, const std::vector<LPPartition> particle_set, const std::vector<double> w, const arma::vec &ybar, const int T, const arma::mat &A_block, const double rho, const double a1, const double a2, const double nu_sigma, const double lambda_sigma, const double eta, const double lambda)
{
  
  if(mi.num_merges > 0){
    LPPartition max_candidate = new Partition(particle_set[current_l]);
    //double max_objective = w[current_l]*total_log_post(max_candidate, a_sigma, nu_sigma) + lambda * Entropy(current_l, max_candidate, particle_set, w);
    double max_objective = w[current_l]*total_log_post(max_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(current_l, max_candidate, particle_set, w);
    LPPartition tmp_candidate = new Partition(particle_set[current_l]); // the running candidate
    double tmp_objective = 0.0;
    
    for(int m_ix = 0; m_ix < mi.num_merges; m_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[current_l]);
      tmp_candidate->Merge(mi.rec_k[m_ix], mi.donor_k[m_ix], ybar, T, A_block, rho, a1, a2, eta);
      //Rcpp::Rcout << "tmp_candidate currently" << endl;
      //tmp_candidate->Print_Partition();
      //tmp_objective = w[current_l]*total_log_post(tmp_candidate, a_sigma, nu_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
      tmp_objective = w[current_l]*total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(current_l, tmp_candidate, particle_set, w);
      if(tmp_objective > max_objective){
        //Rcpp::Rcout << "[best merge]: max_candidate." << endl;
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
      }
    }
    candidate->Copy_Partition(max_candidate);
  } else{
    candidate->Copy_Partition(particle_set[current_l]);
  }
  
}

bool sanity_check(LPPartition partition){
  bool flag = true;
  int K = partition->K;
  int n = partition->nObs;

  int running_count = 0;
  for(int k = 0; k < K; k++){
    running_count += partition->cluster_config[k];
  }
  if(running_count != n){
    flag = false;
  } else{
    flag = true;
  }
  

  
  return flag;
}

