//
//  kmeans.cpp
//  
//
//  Created by Sameer Deshpande on 8/17/19.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"
#include <vector>
#include <ctime>

using namespace arma;
using namespace std;


using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List kmeans_particle(arma::mat Y,
                  const arma::mat A_block,
                  Rcpp::List gamma_init,
                  const int max_split = 5,
                  const double a1 = 1.0,
                  const double a2 = 1.0,
                  const double nu_sigma = 3,
                  const double lambda_sigma = 1,
                  const double rho = 0.99,
                  const double eta = 1.0)
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
  LPPartition gamma_new = new Partition(gamma_0);
  
  // we will run k-means for K = 2, ..., 10
  int split_k = 0; // we will only ever start with gamma_init being the partition with a single cluster
  
  std::vector<LPPartition> particle_set(max_split-1); // holds the results of kmeans after finding connected components
  std::vector<LPPartition> init_particle_set(max_split-1); // hold the actual results of kmeans
  for(int l = 0; l < max_split-1; l++){
    init_particle_set[l] = new Partition(gamma_0);
    particle_set[l] = new Partition(gamma_0);
  }
  
  int n_k = 1;
  
  // initialize stuff for kmeans repeat
  arma::mat U = arma::zeros<mat>(n_k, 1); // holds the data passed to k-means
  arma::mat means = arma::zeros<mat>(n_k,1); // holds data passed to k-means
  double min_score = 0.0; // hold sum-of-squared distances from observations to cluster means
  bool status = true; // tells us whether k-means was successful
  
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // store submatrices of A_block. used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
  
  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters

  
  std::vector<double> time_vec;
  
  time_t tp;
  int time1;
  int time2;
  int l = 0;
  
  
  
  for(int num_splits = 2; num_splits <= max_split; num_splits++){
    //Rcpp::Rcout << "Starting num_splits = " << num_splits << endl;
    n_k = gamma_0->cluster_config[0];
    U.set_size(n_k,1);
    means.set_size(1, num_splits); // one thing to do might be to change to means.zeros(1, num_split).
    min_score = 0.0;
    status = true;
    cluster_dist.set_size(num_splits);
    cluster_dist_indices.set_size(num_splits);
    
    // now we actually do k-means
    for(int ii = 0; ii < n_k; ii++){
      U(ii,0) = ybar(ii); // we are running kmeans on the MLE's
      if(U(ii,0) != U(ii,0)){
        Rcpp::Rcout << "[kmeans]:  possible non-finite value in U" << endl;
      }
    }
    time1 = time(&tp);
    //Rcpp::Rcout << "About to start kmeans repeat" << endl;
    kmeans_repeat(U, means, status, min_score, num_splits, 1000); // run k-means with 5000 random restarts
    if(status == false){
      //Rcpp::Rcout << "kmeans failed!!" << endl;
    } else{
      init_new_clusters.clear();
      init_new_clusters.resize(num_splits);
      for(int i = 0; i < n_k; i++){
        cluster_dist.zeros();
        for(int new_k = 0; new_k < num_splits; new_k++){
          cluster_dist(new_k) = arma::norm(U.row(i).t() - means.col(new_k)); // distance from newly created cluster kk to point i
        }
        cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
        init_new_clusters[cluster_dist_indices(0)].push_back(gamma_0->clusters[split_k][i]);
      } // closes loop over the elements of the original cluster to find new subcluster assignment
      
      //for(int nc_ix = 0; nc_ix < num_splits; nc_ix++){
      //  Rcpp::Rcout << "  new cluster " << nc_ix << " of size : " << init_new_clusters[nc_ix].size() << endl;
      //}
      
      tmp_new_clusters.clear();
      for(int kk = 0; kk < init_new_clusters.size(); kk++){
        connected_components.clear();
        A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // submatrix of A_block corresponds to newly discovered sub-cluster
        new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
        for(int new_k = 0; new_k < connected_components.size(); new_k++){
          tmp_new_clusters.push_back(connected_components[new_k]);
        }
      } // closes loop over the new clusters discovered
      //Rcpp::Rcout << "  With num_splits = " << num_splits << "  actually found " << tmp_new_clusters.size() << "new clusters" << endl;
     
      
      new_clusters.clear();
      new_clusters.resize(tmp_new_clusters.size());
      k_star.clear();
      k_star.resize(tmp_new_clusters.size());
      //Rcpp::Rcout << "  About to create vector new_clusters" << endl;
      for(int nc_ix = 0; nc_ix < new_clusters.size(); nc_ix++){
        for(int ii = 0; ii < tmp_new_clusters[nc_ix].size(); ii++){
          new_clusters[nc_ix].push_back(tmp_new_clusters[nc_ix][ii]);
          k_star[nc_ix] = -1;
        }
      }
      // there seems to be a possible memory leak after Split_Merge
      // a quick workaround would be to use the new Partition constructor
      delete gamma_new;
      gamma_new = new Partition(n, init_new_clusters, ybar, T, A_block, rho, a1, a2, eta);
      init_particle_set[l]->Copy_Partition(gamma_new);
      
      delete gamma_new;
      gamma_new = new Partition(n, new_clusters, ybar, T, A_block, rho, a1, a2, eta);
      particle_set[l]->Copy_Partition(gamma_new);
      l++;
      time2 = time(&tp);
      time_vec.push_back(time2 - time1);
    } // closes else checking that kmeans was successful
  } // closes loop over the number of splits
  
  Rcpp::List init_unik_particles_out;
  format_particle_set(init_particle_set, init_unik_particles_out);
  arma::mat init_alpha_hat_particle = arma::zeros<mat>(n, init_particle_set.size());
  for(int l = 0; l < init_particle_set.size();l++){
    for(int i = 0; i < n; i++){
      init_alpha_hat_particle(i,l) = init_particle_set[l]->alpha_hat[i];
    }
  }
  
  Rcpp::List unik_particles_out;
  std::vector<double> log_like(particle_set.size());
  std::vector<double> log_prior(particle_set.size());
  std::vector<double> log_post(particle_set.size());
  format_particle_set(particle_set, unik_particles_out);
  arma::mat alpha_hat_particle = arma::zeros<mat>(n, particle_set.size());
  for(int l = 0; l < particle_set.size();l++){
    for(int i = 0; i < n; i++){
      alpha_hat_particle(i,l) = particle_set[l]->alpha_hat[i];
    }
    log_like[l] = total_log_like(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
    log_prior[l] = total_log_prior(particle_set[l]);
    log_post[l] = total_log_post(particle_set[l], total_ss, T, nu_sigma, lambda_sigma);
  }
  
  //for(int l = 0; l < particle_set.size(); l++) Rcpp::Rcout << time_vec[l] << endl;
  
  //arma::vec time_out(particle_set.size());
  //for(int l = 0; l < particle_set.size(); l++) time_out(l) = time_vec[l];
  
  
/*
  Rcpp::Rcout << "Finished loop" << endl;
  Rcpp::Rcout << "particle_set_size = " << particle_set.size() << endl;
  
  //for(int l = 0; l < particle_set.size(); l++){
  //  Rcpp::Rcout << l << endl;
  //  particle_set[l]->Print_Partition(total_ss, T, nu_sigma, lambda_sigma);
  //}
  
  
  Rcpp::List unik_particles_out;
  Rcpp::List init_unik_particles_out;
  format_particle_set(particle_set, unik_particles_out);
  format_particle_set(init_particle_set, init_unik_particles_out);
  
  Rcpp::Rcout << "Obtained unik particles" << endl;
  
  arma::mat alpha_hat_particle = arma::zeros<mat>(n, particle_set.size());
  arma::mat init_alpha_hat_particle = arma::zeros<mat>(n, init_particle_set.size());
  for(int l = 0; l < particle_set.size();l++){
    for(int i = 0; i < n; i++){
      alpha_hat_particle(i,l) = particle_set[l]->alpha_hat[i];
      init_alpha_hat_particle(i,l) = init_particle_set[l]->alpha_hat[i];
    }
  }


  */
  Rcpp::List results;
  results["Y"] = Y;
  results["particles"] = unik_particles_out;
  results["alpha"] = alpha_hat_particle;
  results["log_like"] = log_like;
  results["log_prior"] = log_prior;
  results["log_post"] = log_post;
  results["init_particles"] = init_unik_particles_out;
  results["init_alpha_hat_particle"] = init_alpha_hat_particle;
  results["time"] = time_vec;
  return(results);
  
  
  
}
