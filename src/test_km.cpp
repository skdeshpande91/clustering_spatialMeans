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
                   const double split_frac = 0.1)
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
  
  // for testing the new connected components code we need the following
  std::vector<std::vector<int> > test_components;

  
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
        
        status = arma::kmeans(means, U.t(), num_splits, random_subset, 10, false);
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
          
          Rcpp::Rcout << "k-means results: found " << init_new_clusters.size() << " sub-clusters " << endl;
          for(int kk = 0; kk < init_new_clusters.size(); kk++){
            Rcpp::Rcout << "  sub-cluster " << kk << " of size " << init_new_clusters[kk].size()<< " : " ;
            //for(int ii = 0; ii < init_new_clusters[kk].size() ; ii++){
            //  Rcpp::Rcout << " " << init_new_clusters[kk][ii] ;
            //}
            Rcpp::Rcout << endl;
          }
          
          // now loop over init_new_clusters and find the connected components
          tmp_new_clusters.clear();
          for(int kk = 0; kk < init_new_clusters.size(); kk++){
            Rcpp::Rcout << " init_sub_cluster " << kk << " of size " << init_new_clusters[kk].size() ;
            // let's use new_Connected_Components now
            connected_components.clear();
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to the newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            // now connected_components holds all of the connected components of the newly discovered sub-cluster
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
            // now tmp_new_clusters is of length connected_components.size() and each element is a new cluster that we want to work with
            // below is the old code, that gives equivalent answer to new_Connected_Components
            /*
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
            //Rcpp::Rcout << " is broken into " << connected_components.size() << " pieces of size " << endl;
            Rcpp::Rcout << " is broken into " << endl;
            for(int new_k = 0; new_k < connected_components.size() ; new_k++){
              //Rcpp::Rcout << " " << connected_components[new_k].size() ;
              Rcpp::Rcout << "   component " << new_k << " : ";
              for(int ii = 0; ii < connected_components[new_k].size(); ii++){
                Rcpp::Rcout << " " << connected_components[new_k][ii];
              }
              Rcpp::Rcout << endl;
            }
            */
            /*
            // let's try the new_Connected_Components function
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]);
            test_components.clear();
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], test_components);
            //new_Connected_Components(A_tmp,init_new_clusters[kk].size(), init_new_clusters[kk], test_components);
            Rcpp::Rcout << "new connected components are" << endl;
            for(int new_k = 0; new_k < test_components.size(); new_k++){
              Rcpp::Rcout << "component " << new_k << " : " ;
              for(int ii = 0; ii < test_components[new_k].size(); ii++){
                Rcpp::Rcout << " " << test_components[new_k][ii] ;
              }
              Rcpp::Rcout << endl;
            }
           */
          } // closes loop over the number of clusters found by k-means
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
        } // closes if/else checking that k-means succeeded
      } // closes loop over number of potential splits
    } // closes if statement checking that cluster can be split
    
    
    
  } // closes loop over the clusters of gamma_l

  // need to return something
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}
