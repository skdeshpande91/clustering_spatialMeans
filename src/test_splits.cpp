//
//  test_spec_split.cpp
//  
//
//  Created by Sameer Deshpande on 3/13/19.
//

#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "partition.h"
#include "partition_functions.h"
#include "various_functions.h"

using namespace arma;
using namespace std;

using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Rcpp::List test_spec_split(arma::vec ybar, const int T,
                           const arma::mat A_block,
                           Rcpp::List gamma_init,
                           const double a1 = 1.0, const double a2 = 1.0,
                           const double rho = 0.99, const double lambda = 1.0,
                           const double eta = 1.0, const int reps = 1000)
{
  int n = ybar.size();
  LPPartition gamma_l = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  
  split_info si;
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int n_k = 1;
  double min_score = 0.0;
  int split_k = 0;
  
  // initialize all of the stuff for spectral clustering here
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
  bool status = true; // tells us that kmeans was successful
  
  arma::vec cluster_dist = arma::zeros<vec>(1); // holds distance to cluster means
  arma::uvec cluster_dist_indices(1); // used to sort cluster_dist
  
  std::vector<std::vector<int> > init_new_clusters; // stores output of kmeans. NOTE: these clusters may not be connected
  arma::mat A_tmp = arma::zeros<arma::mat>(n_k,n_k);// stores submatrix of A_block, used to determine connected components
  std::vector<std::vector<int> > connected_components; // connected components of individual sub-clusters discovered
  std::vector<std::vector<int> > tmp_new_clusters;// used for figuring out the connected components of new clusters
  std::vector<std::vector<int> > new_clusters; // holds the final new clusters
  std::vector<int> k_star; // actual nearest neighbors of newly formed clusters

  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    split_k = k;
    max_splits = 5;
    if(n_k > 1){
      if(sqrt(n_k) < 5) max_splits = ceil(sqrt(n_k));
      else max_splits = 5;
      for(int num_splits = 2; num_splits <= max_splits; num_splits++){
        
        // re-size everything
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
        means.set_size(num_splits,num_splits);
        status = true;
        min_score = 0.0;
        cluster_dist.set_size(num_splits);
        cluster_dist_indices.set_size(num_splits);
        
        // Spectral Clustering Begins //
        A_block_k = Submatrix(A_block, n_k, n_k, gamma_l->clusters[split_k], gamma_l->clusters[split_k]);
        for(int i = 0; i  < n; i++){
          alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[split_k][i]];
        }
        alpha_dist = Distance_matrix(alpha_hat_cluster,n_k); // distance matrix
        alpha_sim = exp(-1.0 * square(alpha_dist)/(2 * arma::var(alpha_hat_cluster))); // similarity matrix
        W_alpha_cl = I_k + alpha_sim % A_block_k; // ensure non-adjacent indices have 0 similarity
        Dinv_sqrt = arma::diagmat(1/sqrt(arma::sum(W_alpha_cl,1)));
        L = I_k - Dinv_sqrt * W_alpha_cl * Dinv_sqrt;
        arma::eig_sym(eigval, eigvec, L);
        U = eigvec.cols(0,num_splits-1);
        U = arma::diagmat(1/sqrt(arma::sum(arma::square(U),1))) * U;
        // now run kmeans_repeat
        kmeans_repeat(U, means, status, min_score, num_splits, reps);
        
        if(status == false){
          Rcpp::Rcout << "kmeans failed!" << "k = " << k << "num_splits = " << num_splits << endl;
        } else{
          Rcpp::Rcout << "kmeans succeeded!" << endl;
          // do more stuff
          init_new_clusters.clear();
          init_new_clusters.resize(num_splits);
          
          // loop over rows of U to figure out to which cluster it belongs
          for(int i = 0; i < n_k; i++){
            cluster_dist.zeros();
            for(int new_k = 0; new_k < num_splits; new_k++){
              cluster_dist(new_k) = arma::norm(U.row(i).t() - means.col(new_k));
            } // closes loop computing distance of U.row(i) to centroid new_k
            cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
            init_new_clusters[cluster_dist_indices(0)].push_back(gamma_l->clusters[split_k][i]);
          } // closes loop over rows of U
          
          // now loop over init_new_clusters and find the connected components
          tmp_new_clusters.clear();
          for(int kk = 0; kk < init_new_clusters.size(); kk++){
            connected_components.clear();
            A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // submatrix of A_block corresponding to newly discovered sub-cluster
            new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
            for(int new_k = 0; new_k < connected_components.size(); new_k++){
              tmp_new_clusters.push_back(connected_components[new_k]);
            }
          } // closes loop over the number of clusters found by k-means
          // at this point tmp_new_clusters contains all of the new sub-clusters
          
          // get_subcluster_neighbor finds the neighbors of each sub-cluster and arranges them by distance to nearest neighbor
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          // quick printout summarizing the sub-clusters discovered
          Rcpp::Rcout << "  number of subclusters: " << new_clusters.size() << " min score = " << min_score << endl;
          for(int new_k = 0; new_k < new_clusters.size(); new_k++) Rcpp::Rcout << " " << new_clusters[new_k].size();
          Rcpp::Rcout << endl;
          
          // write the information about the split to si
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
        } // closes if/else checking that kmeans_repeat succeeded.
      } // closes loop over number of possible splits
    } // closes if checking that there is at least one element in the cluster
    
    
    
    
  } // closes loop over the clusters
  
  // need to return something
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);

}

// [[Rcpp::export]]
Rcpp::List test_tail_split(arma::vec ybar, const int T,
                           const arma::mat A_block,
                           Rcpp::List gamma_init,
                           const double a1 = 1.0, const double a2 = 1.0,
                           const double rho = 0.99, const double lambda = 1.0,
                           const double eta = 1.0, const double tail_frac = 0.05)
{
  
  int n = ybar.size();
  LPPartition gamma_l = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  
  split_info si;
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int n_k = 1;
  double min_score = 0.0;
  int split_k = 0;

  arma::vec alpha_hat_cluster = arma::zeros<arma::vec>(1);
  arma::uvec alpha_hat_indices(1);
  
  std::vector<int> left_tail;
  std::vector<int> right_tail;
  std::vector<int> center;
  
  std::vector<int> sub_left_tail; // for sub-vectors of left-tail
  std::vector<int> sub_right_tail; // for sub-vectors of right-tail
  std::vector<int> remain; // contains all of the indicies that are not removed from the cluster
  
  std::vector<std::vector<int> > init_new_clusters; // holds the initial sub-clusters
  arma::mat A_tmp; // submatrix corresponding to an initial subcluster
  std::vector<std::vector<int> > connected_components; // connected components of each initial individual sub-cluster
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily holds the new connected sub-clusters
  std::vector<std::vector<int> > new_clusters; // the final sub-clusters arranged in order
  std::vector<int> k_star; // nearest neighbors of newly formed subclusters
  
  for(int k = 0; k < orig_K; k++){
    //Rcpp::Rcout << " k = " << k << endl;
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){
      alpha_hat_cluster.clear();
      alpha_hat_indices.clear();
      alpha_hat_cluster.set_size(n_k); // holds the alpha-hat values in the cluster
      alpha_hat_indices.set_size(n_k);
      for(int i = 0; i < n_k; i++){
        alpha_hat_cluster(i) = gamma_l->alpha_hat[gamma_l->clusters[k][i]];
      }
      alpha_hat_indices = arma::sort_index(alpha_hat_cluster, "ascend");
      left_tail.clear();
      right_tail.clear();
      center.clear();
      for(int i = 0; i < ceil(n_k*tail_frac);i++){
        left_tail.push_back(gamma_l->clusters[k][alpha_hat_indices(i)]);
        right_tail.push_back(gamma_l->clusters[k][alpha_hat_indices(n_k - 1 - i)]);
      }
      for(int i = ceil(n_k*tail_frac); i < n_k - ceil(n_k*tail_frac);i++){
        center.push_back(gamma_l->clusters[k][alpha_hat_indices(i)]);
      }
      
      // now try to remove only the left-tail
      sub_left_tail.clear();
      sub_left_tail.push_back(left_tail[0]);
      
      for(int i = 1; i < left_tail.size(); i++){
        //Rcpp::Rcout << " i = " << i << endl;
        init_new_clusters.clear();
        sub_left_tail.push_back(left_tail[i]); // we are now trying to remove the smallest i elements from cluster
        init_new_clusters.push_back(sub_left_tail); // this subset of left-tail is our first new sub-cluster
        
        // need to figure out what elements are staying in the original cluster
        remain.clear();
        
        for(int ii = i+1; ii < left_tail.size(); ii++) remain.push_back(left_tail[ii]); // all remaining elements of left tail
        for(int ii = 0; ii < center.size(); ii++) remain.push_back(center[ii]); // all elements of center
        for(int ii = 0; ii < right_tail.size(); ii++) remain.push_back(right_tail[ii]); // all elements of right tail
        
        init_new_clusters.push_back(remain); // everything not in sub_left_tail forms our second new sub-cluster
        //Rcpp::Rcout << "    formed the init_new_clusters." << endl;
        // we now need to loop over init_new_clusters and find the connected components
        tmp_new_clusters.clear();
        for(int kk = 0; kk < init_new_clusters.size(); kk++){
          connected_components.clear();
          A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to newly discovered sub-cluster
          new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++) tmp_new_clusters.push_back(connected_components[new_k]);
        } // closes loop over the initial set of sub-clusters
        // at this point, tmp_new_clusters contains all of the new subclusters and they are each connected
        // get_subcluster_neighbor finds the neighbors of each of these subclusters and arranges them by distance to nearest neighbor
        new_clusters.clear();
        new_clusters.resize(tmp_new_clusters.size());
        k_star.clear();
        k_star.resize(tmp_new_clusters.size());
        get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
        
        // quick printout summarizing the sub-clusters discovered
        //Rcpp::Rcout << "  number of sub-clusters: " << new_clusters.size() << endl;
        //for(int new_k = 0; new_k < new_clusters.size(); new_k++) Rcpp::Rcout << " " << new_clusters[new_k].size();
        //Rcpp::Rcout << endl;
        
        si.split_k.push_back(k);
        si.new_clusters.push_back(new_clusters);
        si.nearest_neighbor.push_back(k_star);
        si.num_splits++;
      } // closes loop removing stuff from the left tail
      
      // now remove stuff from the right tail only
      sub_left_tail.clear();
      sub_right_tail.clear();
      sub_right_tail.push_back(right_tail[0]);
      for(int i = 0; i < right_tail.size(); i++){
        init_new_clusters.clear();
        sub_right_tail.push_back(right_tail[i]);
        init_new_clusters.push_back(sub_right_tail);
        
        // need to figure out what elements are staying in the original cluster
        remain.clear();
        
        for(int ii = i+1; ii < right_tail.size(); ii++) remain.push_back(right_tail[ii]); // all remaining elements of right tail
        for(int ii = 0; ii < center.size(); ii++) remain.push_back(center[ii]); // all elements of center
        for(int ii = 0; ii < left_tail.size(); ii++) remain.push_back(left_tail[ii]); // all elements of left tail
        
        init_new_clusters.push_back(remain); // everything not in sub_left_tail forms our second new sub-cluster
        //Rcpp::Rcout << "    formed the init_new_clusters." << endl;
        // we now need to loop over init_new_clusters and find the connected components
        tmp_new_clusters.clear();
        for(int kk = 0; kk < init_new_clusters.size(); kk++){
          connected_components.clear();
          A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to newly discovered sub-cluster
          new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++) tmp_new_clusters.push_back(connected_components[new_k]);
        } // closes loop over the initial set of sub-clusters
        // at this point, tmp_new_clusters contains all of the new subclusters and they are each connected
        // get_subcluster_neighbor finds the neighbors of each of these subclusters and arranges them by distance to nearest neighbor
        new_clusters.clear();
        new_clusters.resize(tmp_new_clusters.size());
        k_star.clear();
        k_star.resize(tmp_new_clusters.size());
        get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
        
        // quick printout summarizing the sub-clusters discovered
        //Rcpp::Rcout << "  number of sub-clusters: " << new_clusters.size() << endl;
        //for(int new_k = 0; new_k < new_clusters.size(); new_k++) Rcpp::Rcout << " " << new_clusters[new_k].size();
        //Rcpp::Rcout << endl;
        si.split_k.push_back(k);
        si.new_clusters.push_back(new_clusters);
        si.nearest_neighbor.push_back(k_star);
        si.num_splits++;
      }
      
      sub_left_tail.clear();
      sub_right_tail.clear();
      
      sub_left_tail.push_back(left_tail[0]);
      sub_right_tail.push_back(right_tail[0]);
      for(int i = 0; i < left_tail.size(); i++){
        init_new_clusters.clear();
        sub_left_tail.push_back(left_tail[i]);
        sub_right_tail.push_back(right_tail[i]);
        init_new_clusters.push_back(sub_left_tail); // first new sub-cluster is from the left
        init_new_clusters.push_back(sub_right_tail); //second new sub-cluster is from the right
        // need to figure out what elements are staying in the original cluster
        remain.clear();
        
        for(int ii = i+1; ii < left_tail.size(); ii++) remain.push_back(left_tail[ii]); // all remaining elements of left tail
        for(int ii = 0; ii < center.size(); ii++) remain.push_back(center[ii]); // all elements of center
        for(int ii = i+1; ii < right_tail.size(); ii++) remain.push_back(right_tail[ii]); // all  remaining elements of right tail
        
        init_new_clusters.push_back(remain); // everything not in sub_left_tail forms our second new sub-cluster
        //Rcpp::Rcout << "    formed the init_new_clusters." << endl;
        // we now need to loop over init_new_clusters and find the connected components
        tmp_new_clusters.clear();
        for(int kk = 0; kk < init_new_clusters.size(); kk++){
          connected_components.clear();
          A_tmp = Submatrix(A_block, init_new_clusters[kk].size(), init_new_clusters[kk].size(), init_new_clusters[kk], init_new_clusters[kk]); // sub-matrix of A_block corresponding to newly discovered sub-cluster
          new_Connected_Components(A_tmp, init_new_clusters[kk].size(), init_new_clusters[kk], connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++) tmp_new_clusters.push_back(connected_components[new_k]);
        } // closes loop over the initial set of sub-clusters
        // at this point, tmp_new_clusters contains all of the new subclusters and they are each connected
        // get_subcluster_neighbor finds the neighbors of each of these subclusters and arranges them by distance to nearest neighbor
        new_clusters.clear();
        new_clusters.resize(tmp_new_clusters.size());
        k_star.clear();
        k_star.resize(tmp_new_clusters.size());
        get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
        
        // quick printout summarizing the sub-clusters discovered
        //Rcpp::Rcout << "  number of sub-clusters: " << new_clusters.size() << endl;
        //for(int new_k = 0; new_k < new_clusters.size(); new_k++) Rcpp::Rcout << " " << new_clusters[new_k].size();
        //Rcpp::Rcout << endl;
        
        si.split_k.push_back(k);
        si.new_clusters.push_back(new_clusters);
        si.nearest_neighbor.push_back(k_star);
        si.num_splits++;
      }
      
      // now remove stuff from both tails simultaneously
    } // closes if/else checking that we can split cluster (n_k > 1)
  }
  
  // need to return something
  List results;
  results["ybar"] = ybar;
  return(results);
}

// [[Rcpp::export]]
Rcpp::List test_border(arma::vec ybar, const int T,
                           const arma::mat A_block,
                           Rcpp::List gamma_init,
                           const double a1 = 1.0, const double a2 = 1.0,
                           const double rho = 0.99, const double lambda = 1.0,
                           const double eta = 1.0)
{
  int n = ybar.size();
  LPPartition gamma_l = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  
  split_info si;
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int n_k = 1;
  double min_score = 0.0;
  int split_k = 0;
  
  arma::mat A_sub = arma::zeros<arma::mat>(1,1); // sub-matrix of A_sub used to see if clusters k and kk are adjacenct
  std::vector<int> adj_cluster(1); // temporarily holds the labels of adjacent clusters
  arma::vec adj_cluster_dist = arma::zeros<arma::vec>(1); // holds distance between cluster k and its adjacent clusters
  arma::uvec adj_cluster_indices(1); // for sorting the adjacent clusters
  
  
  std::vector<int> border(1); // holds indices of the border elements that move from split_k to k
  std::vector<int> remain(1); // holds indices that remain in cluster split_k (i.e. they are not adjacent to block-group k)
  arma::mat A_tmp = arma::zeros<arma::mat>(1,1); // used to find connected components of remain
  std::vector<std::vector<int> > connected_components; // temporariy hold connected components of each sub-cluster
  std::vector<std::vector<int> > new_clusters; // holds the new connected sub-clusters created
  std::vector<int> k_star; // holds the labels of nearest neighbors of new sub-clusters

  if(orig_K > 1){
    for(int k = 0; k < orig_K; k++){
      adj_cluster.clear();
      for(int kk = 0; kk < orig_K; kk++){
        if(kk != k){
          A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[kk], gamma_l->clusters[k], gamma_l->clusters[kk]);
          if(any(vectorise(A_sub) == 1.0)) adj_cluster.push_back(kk); // cluster k and kk are adjacent
        } // closes if checking that kk != k
      } // closes loop over all clusters
      if(adj_cluster.size() > 0){
        adj_cluster_dist.resize(adj_cluster.size());
        adj_cluster_indices.resize(adj_cluster.size());
        for(int kk = 0; kk < adj_cluster.size(); kk++){
          adj_cluster_dist(kk) = abs(gamma_l->alpha_bar[k] - gamma_l->alpha_bar[adj_cluster[kk]]); // distance between cluster k and kk
        }
        adj_cluster_indices = arma::sort_index(adj_cluster_dist, "ascend");
        split_k = adj_cluster[adj_cluster_indices(0)];
        Rcpp::Rcout << "k = " << k << " : will move border elements from " << split_k << " into " << k << endl;
        A_sub = Submatrix(A_block, gamma_l->cluster_config[k], gamma_l->cluster_config[split_k], gamma_l->clusters[k], gamma_l->clusters[split_k]); // submatrix of A_block corresponding to clusters k and its neighbor split_k
        
        border.clear();
        remain.clear();

        for(int i = 0; i < gamma_l->cluster_config[split_k]; i++){
          // anything in A_sub.col(i) is equal to 1 then element i in cluster split_k is a border
          if(any(vectorise(A_sub.col(i)) == 1.0)) border.push_back(gamma_l->clusters[split_k][i]);
          else remain.push_back(gamma_l->clusters[split_k][i]);
        }
        if(remain.size() > 0 & border.size() > 0){ // condition may be a bit redundant
          // first find the connected components of the elements in remain
          new_clusters.clear();
          k_star.clear();
          connected_components.clear();
          A_tmp = Submatrix(A_block, remain.size(), remain.size(), remain, remain);
          new_Connected_Components(A_tmp, remain.size(), remain, connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++){
            new_clusters.push_back(connected_components[new_k]);
            k_star.push_back(-1);
          }
          // up to now new_clusters and k_star don't include the border elements which are moving to cluster k
          new_clusters.push_back(border);
          k_star.push_back(k);
 
          si.split_k.push_back(split_k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
          
          Rcpp::Rcout << "Removing " << border.size() << " indices from cluster " << split_k << " and moving to cluster " << k << endl;
          for(int i = 0; i < border.size(); i++) Rcpp::Rcout << " " << border[i] ;
          Rcpp::Rcout << endl;
          
          
        } // closes if checking that not all elements in split_k are border (if they are don't do anything!)
      } // closes if checking that there are valid border moves
    } // closes loop over the clusters k
  } // check that there are at least two clusters
  
  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
  
}

// [[Rcpp::export]]
Rcpp::List test_island(arma::vec ybar, const int T,
                       const arma::mat A_block,
                       Rcpp::List gamma_init,
                       const double a1 = 1.0, const double a2 = 1.0,
                       const double rho = 0.99, const double lambda = 1.0,
                       const double eta = 1.0, const double island_frac = 0.05)
{
  int n = ybar.size();
  LPPartition gamma_l = new Partition(n, gamma_init, ybar, T, A_block, rho, a1, a2, eta);
  
  split_info si;
  int orig_K = gamma_l->K;
  int max_splits = 5;
  int n_k = 1;
  double min_score = 0.0;
  int split_k = 0;
  double alpha_bar = 0.0; // holds the cluster mean
  arma::vec distance = arma::zeros<arma::vec>(1); // holds distance of each alpha-hat from the overall cluster mean
  arma::uvec indices(1); // used to sort the elements based on how far they are from the cluster mean
  int num_islands = 0; // number of islands we try to create
  
  std::vector<int> island(1);
  std::vector<int> remain(1); // holds indices that remain in cluster split_k (i.e. they are not adjacent to block-group k)
  arma::mat A_tmp = arma::zeros<arma::mat>(1,1); // used to find connected components of remain
  std::vector<std::vector<int> > connected_components; // temporariy hold connected components of each sub-cluster
  std::vector<std::vector<int> > tmp_new_clusters; // will contains connected components of remain and the island
  std::vector<std::vector<int> > new_clusters; // holds the new connected sub-clusters created
  std::vector<int> k_star; // holds the labels of nearest neighbors of new sub-clusters
  
  
  for(int k = 0; k < orig_K; k++){
    n_k = gamma_l->cluster_config[k];
    if(n_k > 1){
      distance.reset();
      distance.set_size(n_k);
      indices.reset();
      indices.set_size(n_k);
      
      alpha_bar = gamma_l->alpha_bar[k];
      for(int i = 0; i < n_k; i++) distance(i) = abs(alpha_bar - gamma_l->alpha_hat[gamma_l->clusters[k][i]]);
      indices = arma::sort_index(distance, "descend");
      num_islands = ceil( (double) n_k * island_frac);
      for(int ix = 0; ix < num_islands; ix++){
        remain.clear();
        island.clear();
        for(int i = 0; i < n_k; i++){
          if(i != ix) remain.push_back(gamma_l->clusters[k][indices(i)]); // includes all but the ix-th furtherst point from cluster mean
          else island.push_back(gamma_l->clusters[k][indices(i)]);
        } // closes loop populating remain and island
        
        // find the connected components of remain
        if(remain.size() > 0 & island.size() > 0){ // this is probably superfluous but it's fine for now
          connected_components.clear();
          tmp_new_clusters.clear();
          A_tmp = Submatrix(A_block, remain.size(), remain.size(), remain, remain);
          new_Connected_Components(A_tmp, remain.size(), remain, connected_components);
          for(int new_k = 0; new_k < connected_components.size(); new_k++){
            tmp_new_clusters.push_back(connected_components[new_k]);
          }
          tmp_new_clusters.push_back(island); // add the island to the mix
          new_clusters.clear();
          new_clusters.resize(tmp_new_clusters.size());
          k_star.clear();
          k_star.resize(tmp_new_clusters.size());
          get_subcluster_neighbor(tmp_new_clusters, new_clusters, k_star, k, gamma_l, T, A_block, rho, a1, a2);
          
          si.split_k.push_back(k);
          si.new_clusters.push_back(new_clusters);
          si.nearest_neighbor.push_back(k_star);
          si.num_splits++;
          
          // quick printout to check progress
          Rcpp::Rcout << "  creating " << new_clusters.size() << " new clusters" << endl;
          for(int new_k = 0; new_k < new_clusters.size(); new_k++){
            Rcpp::Rcout << "     sub cluster " << new_k << " of size " << new_clusters[new_k].size() << " and neighbor " << k_star[new_k] << endl;
          }
        } // closes if/else checking that at least 1 elements remains in k.

      } // closes loop over the possible islands
    } // closes if/else checking that there is more than 1 element in the cluster
  } // closes loop over all clusters


  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}




// [[Rcpp::export]]
Rcpp::List test_moves(arma::vec ybar, const int T,
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
  
  LPPartition max_candidate = new Partition(gamma_0);
  LPPartition tmp_candidate = new Partition(gamma_0);
  
  double max_objective = 0.0;
  double tmp_objective = 0.0;
  
  split_info spec_si;
  split_info km_si;
  split_info tail_si;
  split_info isl_i; // holds information for island moves
  split_info bi; // holds information for border moves
  merge_info mi;
  
  for(int l = 0; l < L; l++){
    Rcpp::Rcout << "Starting particle " << l << endl;
    Rcpp::checkUserInterrupt();
    delete max_candidate;
    max_candidate = new Partition(particle_set[l]);
    max_objective = w[l]*total_log_post(max_candidate, nu_sigma, lambda_sigma) + lambda*Entropy(l, max_candidate, particle_set, w);
    
    get_island(isl_i, particle_set[l], T, A_block, rho, a1, a2, 0.05);
    Rcpp::Rcout << "Got " << isl_i.num_splits << " possible island splits" << endl;
    for(int split_ix = 0; split_ix < isl_i.num_splits; split_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Split_Merge(isl_i.split_k[split_ix], isl_i.new_clusters[split_ix], isl_i.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[l] * total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, tmp_candidate, particle_set,w);
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        Rcpp::Rcout << "updated with island proposal" << endl;
        max_candidate->Print_Partition(nu_sigma, lambda_sigma);
      }
    }
    
    get_border(bi, particle_set[l], T, A_block, rho, a1, a2);
    Rcpp::Rcout << "Got " << bi.num_splits << " possible border moves" << endl;
    for(int split_ix = 0; split_ix < bi.num_splits; split_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Split_Merge(bi.split_k[split_ix], bi.new_clusters[split_ix], bi.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[l] * total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, tmp_candidate, particle_set,w);
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        Rcpp::Rcout << "updated with border move" << endl;
        max_candidate->Print_Partition(nu_sigma, lambda_sigma);
      }
    }
    
    get_merge(mi, particle_set[l], A_block);
    Rcpp::Rcout << "Got " << mi.num_merges << " possible merge moves" << endl;
    for(int m_ix = 0; m_ix < mi.num_merges; m_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Merge(mi.rec_k[m_ix], mi.donor_k[m_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[l]*total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, tmp_candidate, particle_set, w);
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        Rcpp::Rcout << "Updated with merge move" << endl;
        max_candidate->Print_Partition(nu_sigma, lambda_sigma);
      }
    }
    
    get_tail_split(tail_si, particle_set[l], T, A_block, rho, a1, a2, 0.025);
    Rcpp::Rcout << "Got " << tail_si.num_splits << " possible tail splits" << endl;
    for(int split_ix = 0; split_ix < tail_si.num_splits; split_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Split_Merge(tail_si.split_k[split_ix], tail_si.new_clusters[split_ix], tail_si.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
      //tmp_candidate->Print_Partition(nu_sigma, lambda_sigma);
    }

    // test spectral split
    get_spectral_split(spec_si, particle_set[l], T, A_block, rho, a1, a2, 10);
    Rcpp::Rcout << "Got " << spec_si.num_splits << "  possible spectral splits" << endl;
    for(int split_ix = 0; split_ix < spec_si.num_splits; split_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Split_Merge(spec_si.split_k[split_ix], spec_si.new_clusters[split_ix], spec_si.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[l] * total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, tmp_candidate, particle_set,w);
      //tmp_candidate->Print_Partition(nu_sigma, lambda_sigma);
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        Rcpp::Rcout << "Updated the particle with a spectral split!" << endl;
      }
    }

   // now try km splits
    get_km_split(km_si, particle_set[l], T, A_block, rho, a1, a2, 5000);
    for(int split_ix = 0; split_ix < km_si.num_splits; split_ix++){
      delete tmp_candidate;
      tmp_candidate = new Partition(particle_set[l]);
      tmp_candidate->Split_Merge(km_si.split_k[split_ix], km_si.new_clusters[split_ix], km_si.nearest_neighbor[split_ix], ybar, T, A_block, rho, a1, a2, eta);
      tmp_objective = w[l] * total_log_post(tmp_candidate, nu_sigma, lambda_sigma) + lambda * Entropy(l, tmp_candidate, particle_set, w);
      if(tmp_objective > max_objective){
        delete max_candidate;
        max_candidate = new Partition(tmp_candidate);
        max_objective = tmp_objective;
        Rcpp::Rcout << "Updated the particle with a km split!" << endl;
      }
    }
  }

  Rcpp::List results;
  results["ybar"] = ybar;
  return(results);
}
