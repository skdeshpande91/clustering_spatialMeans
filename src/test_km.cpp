//
//  test_km.cpp
//  
//
//  Created by Sameer Deshpande on 2/22/19.
//
#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include "various_functions.h"
#include <vector>

using namespace arma;
using namespace std;


using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List test_km(arma::vec ybar, const int T,
                   const arma::mat A_block,
                   Rcpp::List gamma_init,
                   const double a1 = 1.0, const double a2 = 1.0,
                   const double nu_sigma = 1, const double lambda_sigma = 1,
                   const double rho = 0.99, const double lambda = 1.0,
                   const double eta = 1.0,
                   const int reps = 1000)
{
  int n = ybar.size();
  LPPartition gamma_l = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  
  // eventually we should create a particle set with L = 10
  // look at what splits are proposed and which are accepted
  
  split_info si;
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int split_k = 0;
  int n_k = 1;
  double min_score = 0.0;

  // initialize the stuff needed for k-means clustering
  arma::mat U = arma::zeros<mat>(n_k, 1); // holds the data passed to k-means
  arma::mat means = arma::zeros<mat>(n_k,1); // holds data passed to k-means
  bool status = true; // tells us whether k-means was successful
  
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // stores submatrices of A_block. used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
  
  
  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
  
  for(int k = 0; k < orig_K; k++){
    //Rcpp::Rcout << "[get_km_split]: k = " << k << endl;
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    max_splits = 5;
    if(n_k > 1){
      //if(sqrt(n_k) < 5) max_splits = sqrt(n_k);
      //else max_splits = 5;
      //if(2 >= ceil(split_frac * sqrt(n_k))) max_splits = 3;
      // attempt only up to sqrt(n_k) * split_frac splits.
      //else max_splits = ceil(sqrt(n_k) * split_frac);
      //Rcpp::Rcout << "max_splits = " << max_splits << endl;
      for(int num_splits = 2; num_splits <= max_splits; num_splits++){
        Rcpp::Rcout << "Starting num_splits = " << num_splits;
        U.set_size(n_k, 1);
        means.set_size(1, num_splits);
        status = true;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        //Rcpp::Rcout << "    set U, means, status " << endl;
        
        // do k-means
        for(int ii = 0; ii < n_k; ii++){
          U(ii,0) = gamma_l->alpha_hat[gamma_l->clusters[split_k][ii]];
          if(U(ii,0) != U(ii,0)){
            Rcpp::Rcout << "  possible non-finite value in U" << endl;
            Rcpp::Rcout << "split_k = " << split_k << "  ii = " << ii << endl;
            Rcpp::Rcout << gamma_l->clusters[split_k][ii] << endl;
          }
        }
        //Rcpp::Rcout << "Ready to start k-means" << endl;
        kmeans_repeat(U, means, status, min_score, num_splits, reps);
        
        //status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
        //status = arma::kmeans(means, U.t(), num_splits, random_spread, 10, false);
        
        if(status == false){
          Rcpp::Rcout << "kmeans failed!!" << "n_k = " << n_k << " num_splits = " << num_splits << endl;
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
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to the newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the number of clusters found by k-means
          // now tmp_new_clusters contains all of the new subclusters discovered by kmeans.

          // get_subcluster_neighbor finds the neighbors of each sub-cluster and arranges them by distance to nearest neighbor
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          
          // quick printout summarizing the sub-clusters discovered
          Rcpp::Rcout << "  number of sub-clusters: " << new_clusters.size() << "  min score = " << min_score << endl;
          for(int new_k = 0; new_k < new_clusters.size(); new_k++)  Rcpp::Rcout << " " << new_clusters[new_k].size();
          Rcpp::Rcout << endl;
          
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;

        } // closes if/else checking that k-means succeeded
      } // closes loop over number of potential splits
    } // closes if statement checking that cluster can be split
  } // closes loop over the clusters of gamma_l

  
  
  // now let's run get_km_split. This should yield the same results
  split_info new_si;
  get_km_split(new_si, gamma_l, T, A_block, rho, a1, a2, 500);
  
  
  
  // need to return something
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}

// [[Rcpp::export]]
Rcpp::List test_kmpp(arma::vec ybar, const int T,
                   const arma::mat A_block,
                   Rcpp::List gamma_init,
                   const double a1 = 1.0, const double a2 = 1.0,
                   const double nu_sigma = 1, const double lambda_sigma = 1,
                   const double rho = 0.99, const double lambda = 1.0,
                   const double eta = 1.0,
                   const int reps = 1000)
{
  
  int n = ybar.size();
  LPPartition gamma_l = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  
  // eventually we should create a particle set with L = 10
  // look at what splits are proposed and which are accepted
  
  split_info si;
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int split_k = 0;
  int n_k = 1;
  double min_score = 0.0;
  
  // initialize the stuff needed for k-means clustering
  arma::mat U = arma::zeros<mat>(n_k, 1); // holds the data passed to k-means
  arma::mat means = arma::zeros<mat>(n_k,1); // holds data passed to k-means
  bool status = true; // tells us whether k-means was successful
  
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of k-means. NOTE: these clusters may not be connected.
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k); // stores submatrices of A_block. used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of the individual sub clusters found by k-means
  
  
  std::vector<std::vector<int> > tmp_new_clusters; // used for figuring out the connected components of the new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters
  
  for(int k = 0; k < orig_K; k++){
    //Rcpp::Rcout << "[get_km_split]: k = " << k << endl;
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    max_splits = 5;
    if(n_k > 1){
      //if(sqrt(n_k) < 5) max_splits = sqrt(n_k);
      //else max_splits = 5;
      //if(2 >= ceil(split_frac * sqrt(n_k))) max_splits = 3;
      // attempt only up to sqrt(n_k) * split_frac splits.
      //else max_splits = ceil(sqrt(n_k) * split_frac);
      //Rcpp::Rcout << "max_splits = " << max_splits << endl;
      for(int num_splits = 2; num_splits <= max_splits; num_splits++){
        Rcpp::Rcout << "Starting num_splits = " << num_splits << endl;
        U.set_size(n_k, 1);
        means.set_size(1, num_splits);
        status = true;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        //Rcpp::Rcout << "    set U, means, status " << endl;
        
        // do k-means
        for(int ii = 0; ii < n_k; ii++){
          U(ii,0) = gamma_l->alpha_hat[gamma_l->clusters[split_k][ii]];
          if(U(ii,0) != U(ii,0)){
            Rcpp::Rcout << "  possible non-finite value in U" << endl;
            Rcpp::Rcout << "split_k = " << split_k << "  ii = " << ii << endl;
            Rcpp::Rcout << gamma_l->clusters[split_k][ii] << endl;
          }
        }
        //Rcpp::Rcout << "Ready to start k-means" << endl;
        //kmeans_repeat(U, means, status, min_score, num_splits, reps);
        kmeans_plus_plus(U, means, status, min_score, num_splits, reps);
        
        //status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
        //status = arma::kmeans(means, U.t(), num_splits, random_spread, 10, false);

        if(status == false){
          Rcpp::Rcout << "kmeans failed!!" << "n_k = " << n_k << " num_splits = " << num_splits << endl;
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
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to the newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the number of clusters found by k-means
          // now tmp_new_clusters contains all of the new subclusters discovered by kmeans.
          
          // get_subcluster_neighbor finds the neighbors of each sub-cluster and arranges them by distance to nearest neighbor
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          
          // quick printout summarizing the sub-clusters discovered
          Rcpp::Rcout << "  number of sub-clusters: " << new_clusters.size() << "  min score = " << min_score << endl;
          for(int new_k = 0; new_k < new_clusters.size(); new_k++)  Rcpp::Rcout << " " << new_clusters[new_k].size();
          Rcpp::Rcout << endl;
          

          
          
        } // closes if/else checking that k-means succeeded
      } // closes loop over number of potential splits
    } // closes if statement checking that cluster can be split
  } // closes loop over the clusters of gamma_l

  // need to return something
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}


// [[Rcpp::export]]
Rcpp::List test_km_particle(arma::vec ybar, const int T,
                     const arma::mat A_block,
                     const int L,
                     Rcpp::List gamma_init,
                     const double a1 = 1.0, const double a2 = 1.0,
                     const double nu_sigma = 1, const double lambda_sigma = 1,
                     const double rho = 0.99, const double lambda = 1.0,
                     const double eta = 1.0,
                     const int reps = 1000)
{
  
  int n = ybar.size();
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  // create a particle set
  std::vector<LPPartition> particle_set(L);
  std::vector<double> w(L);
  for(int l = 0; l < L; l++){
    particle_set[l] = new Partition(gamma_0);
    w[l] = 1.0/( (double) L);
  }
  // eventually we should create a particle set with L = 10
  // look at what splits are proposed and which are accepted
  
  LPPartition tmp_candidate = new Partition(gamma_0);
  LPPartition max_candidate = new Partition(gamma_0);
  double max_objective = 0.0;
  double tmp_objective = 0.0;
  split_info km_si;
  
  int num_new_clusters = 0;
  
  
  
  
  for(int l = 0; l < L; l++){
    Rcpp::Rcout << "Starting particle " << l << endl;
    Rcpp::checkUserInterrupt();
    delete max_candidate;
    max_candidate = new Partition(particle_set[l]);
    max_objective = w[l]*total_log_post(max_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, max_candidate, particle_set, w);

    get_km_split(km_si, particle_set[l], T, A_block, rho, a1, a2, 5000);
    
    // now loop through the elements in km_si, evaluate the candidates and print out their objective values.
    // don't immediately call best_split yet though
    Rcpp::Rcout << "Particle " << l << " is :" << endl;
    //particle_set[l]->Print_Partition(nu_sigma, lambda_sigma);
    
    Rcpp::Rcout << "Got " << km_si.num_splits << " possible splits" << endl;
    for(int split_ix = 0; split_ix < km_si.num_splits; split_ix++){
      //Rcpp::Rcout << "Split candidate " << split_ix << endl;
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Split_Merge(km_si.split_k[split_ix], km_si.new_clusters[split_ix], km_si.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[l]*total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, tmp_candidate, particle_set, w);
      Rcpp::Rcout << "tmp_obj = " << tmp_objective << "    max_obj = " << max_objective << endl;
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        Rcpp::Rcout << "Updated the particle!" << endl;
      }
    }
    max_candidate->Print_Partition(nu_sigma, lambda_sigma);
  }
  
  
  
  
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}
