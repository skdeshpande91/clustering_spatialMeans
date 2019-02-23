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
  
  split_info si;
  
  
  
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
      if(sqrt(n_k) < 5) max_splits = sqrt(n_k);
      else max_splits = 5;
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
        
        
        kmeans_repeat(U, means, status, num_splits, reps);
        
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
            //Rcpp::Rcout << " initially sub-cluster " << kk << " has size " << init_new_clusters[kk].size() ;
            // let's use new_Connected_Components now
            connected_components.clear();
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to the newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            //Rcpp::Rcout << "    has " << connected_components.size() << "  connected components." << endl;
            // now connected_components holds all of the connected components of the newly discovered sub-cluster
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the number of clusters found by k-means
          // now tmp_new_clusters contains all of the new subclusters discovered by kmeans.

          // Now we can arrange them
          
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          Rcpp::Rcout << "  number of sub-clusters: " << new_clusters.size() << endl;
          for(int new_k = 0; new_k < new_clusters.size(); new_k++){
            //Rcpp::Rcout << "  subcluster " << new_k << " of size " << new_clusters[new_k].size() << " and neighbor " << k_star[new_k] << endl;
            Rcpp::Rcout << " " << new_clusters[new_k].size();
          }
          Rcpp::Rcout << endl;
          
          //Rcpp::Rcout << "   got connected components following k-means split " << endl;
/*
          for(int new_k = 0; new_k < tmp_new_clusters.size(); new_k++){
            Rcpp::Rcout << "new cluster " << new_k << " : " << endl;
            for(int ii = 0; ii < tmp_new_clusters[new_k].size(); ii++){
              Rcpp::Rcout << " " << tmp_new_clusters[new_k][ii] ;
            }
            Rcpp::Rcout << endl;
          }
*/
        } // closes if/else checking that k-means succeeded
      } // closes loop over number of potential splits
    } // closes if statement checking that cluster can be split
    
    
    
  } // closes loop over the clusters of gamma_l

  // need to return something
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}
