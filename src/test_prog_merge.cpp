//
//  test_prog_merge.cpp
//  Function that is meant to test the progressive merges
//
//  Created by Sameer Deshpande on 2/2/19.
//

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <vector>
#include "partition.h"
#include "various_functions.h"
#include "partition_functions.h"


using namespace arma;
using namespace std;

using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List test_prog_merge(arma::vec ybar,
                           const int T,
                           const arma::mat A_block,
                           const int L,
                           Rcpp::List gamma_init,
                           const double a1 = 1.0,
                           const double a2 = 1.0,
                           const double nu_sigma = 1,
                           const double lambda_sigma = 1,
                           const double rho = 0.99,
                           const double lambda = 1.0,
                           const double eta = 1.0,
                           const int max_iter = 10,
                           const double eps = 1e-3,
                           const double split_frac = 0.1)
{
  int n = ybar.size();
  Rcpp::Rcout << "n = " << endl;
  LPPartition gamma_0 = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  Rcpp::Rcout << "Created gamma_0" << endl;
  
  //gamma_0->Print_Partition(nu_sigma, lambda_sigma);
  
  // get some splits
  split_info si; // holds information for some splits
  
  get_km_split(si, gamma_0, T, A_block, rho, a1, a2, split_frac);
  /*
  Rcpp::Rcout << "Proposed " << si.num_splits << " potential splits" << endl;
  for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
    Rcpp::Rcout << "  Split cluster " << si.split_k[split_ix] << " into " << si.new_clusters[split_ix].size() << " new clusters" << endl;
    for(int nc_ix = 0; nc_ix < si.new_clusters[split_ix].size(); nc_ix++){
      Rcpp::Rcout << "  new cluster " << nc_ix << "  of size " << si.new_clusters[split_ix][nc_ix].size() << " and nearest neighbor " << si.nearest_neighbor[split_ix][nc_ix] << ":" << endl;
      for(int i = 0; i < si.new_clusters[split_ix][nc_ix].size(); i++){
        Rcpp::Rcout << si.new_clusters[split_ix][nc_ix][i] << " ";
      } // closes loop over sub-cluster members
      Rcpp::Rcout << endl;
    } // closes loop over new sub-clusters
  } // closes loop over the proposed splits
  */
  
  LPPartition split_candidate = new Partition(gamma_0);
  
  // now try to find the best split candidates
  // instead of calling the best_split function, let's just write out the code
  
  std::vector<int> k_star; // holds the indices of nearest neighbors
  int num_new_clusters = 0;
  int split_k = 0;

  bool sanity_flag = true;
  //for(int split_ix = 0; split_ix < si.num_splits; split_ix++){
  for(int split_ix = 0; split_ix < 1; split_ix++){
    Rcpp::Rcout << "Starting to consider split " << split_ix << endl;
    
    split_k = si.split_k[split_ix];
    num_new_clusters = si.new_clusters[split_ix].size();
    k_star.clear();
    k_star.resize(num_new_clusters, -1);
    //Rcpp::Rcout << "Cluster " << split_k << " split into " << num_new_clusters << " sub-clusters" << endl;
    //Rcpp::Rcout << "  Nearest neighbors are:";
    //for(int nc_ix = 0; nc_ix < num_new_clusters; nc_ix++){
    //  Rcpp::Rcout << " " << si.nearest_neighbor[split_ix][nc_ix] ;
    //}
    //Rcpp::Rcout << endl;
    //Rcpp::Rcout << " Changing k_star now" << endl;
    
    for(int nc_ix = 0; nc_ix < num_new_clusters;nc_ix++){
      if(si.nearest_neighbor[split_ix][nc_ix] != -1){
        k_star[nc_ix] = si.nearest_neighbor[split_ix][nc_ix]; // as loop progresses, we change one nearest neighbor at a time.
        delete split_candidate;
        split_candidate = new Partition(gamma_0);
        split_candidate->Split_Merge(split_k, si.new_clusters[split_ix], k_star, ybar, T, A_block, rho, a1, a2, eta);
        split_candidate->Print_Partition(nu_sigma, lambda_sigma);
      } // closes if checking that merging a new sub-cluster
      
      
      //Rcpp::Rcout << " sub-cluster " << nc_ix << " k_star is now :";
      //for(int nn = 0; nn < num_new_clusters; nn++){
      //  Rcpp::Rcout << " " << k_star[nn];
      //}
      //Rcpp::Rcout << endl;
      
    }
    
  } // closes loop over the potential splits
  
  
  
  
  
  // This is a placeholder for now.
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}

