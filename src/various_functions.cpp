//
//  various_functions.cpp
//  
//
//  Created by Sameer Deshpande on 10/29/18.
//

#include "various_functions.h"


arma::mat Submatrix(arma::mat M, int n_row, int n_col, std::vector<int> row_index, std::vector<int> col_index)
{
  arma::mat N(n_row, n_col);
  //Rcpp::Rcout << "[Submatrix]: Made N with " << n_row << " rows and " << n_col << "columns " << arma::endl;
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_col; j++){
      //Rcpp::Rcout << row_index[i] << "," << col_index[j] << arma::endl;
      N(i,j) = M(row_index[i], col_index[j]);
    }
  }
  return N;
}
// will need to update these potentially

void Connected_Components(const arma::mat M, const int n, int* components, int* count)
{
  // Mark all the vertices as not visited
  bool *visited = new bool[n];
  *count = 0;
  for(int v = 0; v < n; v++)
    visited[v] = false;
  //Rcpp::Rcout << "[Connected_Components]: Starting main loop" << arma::endl;
  for(int v = 0; v < n; v++)
  {
    //Rcpp::Rcout << "visiting " << v << arma::endl;
    if(visited[v] == false)
    {
      DFSUtil(M, n, v, visited, components, count);
      (*count)++;
    }
  }
  delete[] visited;
}

// to be used in the split functions.
// reads in the sub-cluster found by the split (init_components) and forms the connected components (in the vector of vectors components)
void new_Connected_Components(const arma::mat &M, const int n, std::vector<int> &init_components, std::vector<std::vector<int> >&components)
{
  //Rcpp::Rcout << "[new_Connected_Components]: starting" << std::endl;
  int* count = new int;
  int* tmp_components = new int[n];
  bool* visited = new bool[n];
  
  *count = 0;
  for(int v = 0; v < n; v++){
    visited[v] = false;
  }
  //Rcpp::Rcout << "[new_Connected_Components] : about to start DFS" << std::endl;
  for(int v = 0; v < n; v++){
    if(visited[v] == false){
      DFSUtil(M,n,v,visited, tmp_components,count);
      (*count)++;
    }
  }
  //Rcpp::Rcout << "[new_Connected_Components]: Finished DFS" << std::endl;
  
  // use tmp_components to indx init_components
  components.clear();
  components.resize(*count);
  //int cluster_id = 0;
  for(int i = 0; i < n; i++){
    components[tmp_components[i]].push_back(init_components[i]);
  }
  delete count;
  delete[] tmp_components;
  delete[] visited;
}



/*
void new_Connected_Components(arma::mat M, int n, std::vector<int> components, int count){
  // re-size components
  components.resize(n,0);
  count = 0;
  std::vector<bool> visited(n, false);
  for(int v = 0; v < n; v++){
    if(visited[v] == false){
      new_DFSUtil(M, n, v, visited, components, count);
      count++;
    }
  }
}
*/



void DFSUtil(const arma::mat M, const int n, int v, bool* visited, int* components, int* count)
{
  visited[v] = true;
  components[v] = *count;
  for(int i = 0; i < n; i++)
    if(M(v,i) == 1.0)
      if(visited[i] == false)
        DFSUtil(M, n, i, visited, components, count);
}
/*
void new_DFSUtil(arma::mat M, int n, int v, std::vector<bool> visited, std::vector<int> components, int count)
{
  visited[v] = true;
  components[v] = count;
  for(int i = 0; i < n; i++){
    if(M(v,i) == 1.0){
      if(visited[i] == false){
        new_DFSUtil(M, n, i, visited, components, count);
      }
    }
  }
}
*/
arma::mat Distance_matrix(arma::vec alpha_hat, int nObs){
  arma::mat dist = arma::zeros<arma::mat>(nObs, nObs);
  for(int i = 0; i < nObs; i++){
    for(int j = i+1; j < nObs; j++){
      dist(i,j) = abs(alpha_hat[i] - alpha_hat[j]);
      dist(j,i) = dist(i,j);
    }
  }
  return dist;
}







std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain, const arma::mat &A_block){
  std::vector<std::vector<int> > remain_clusters; 
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily used to build the connected components.
  remain_clusters.push_back(std::vector<int>(1,remain[0]));
  arma::mat A_tmp; 
  for(int ii = 1; ii < remain.size(); ii++){
    // consider the next element in remain and check if it's connected to any of the remain_clusters
    tmp_new_clusters.push_back(std::vector<int>(1, remain[ii]));
    // loop over the existing connected components
    for(int conncomp = 0; conncomp < remain_clusters.size(); conncomp++){
      A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), remain_clusters[conncomp].size(), tmp_new_clusters[0], remain_clusters[conncomp]);
      if(any(vectorise(A_tmp) == 1.0)){ // something in remain_clusters[conncomp] is adjacent to remain[ii] so we should combine these clusters
        for(int ix = 0; ix < remain_clusters[conncomp].size(); ix++){
          tmp_new_clusters[0].push_back(remain_clusters[conncomp][ix]);
        }
      } else{ // remain_clusters[conncomp] remains its own distinct component
        tmp_new_clusters.push_back(remain_clusters[conncomp]);
      }
      remain_clusters[conncomp].clear();
    } // closes loop over elements of remain_clusters
    // update remain_clusters: copy tmp_new_clusters
    remain_clusters.clear();
    for(int cc = 0; cc < tmp_new_clusters.size(); cc++){
      remain_clusters.push_back(tmp_new_clusters[cc]);
      tmp_new_clusters[cc].clear();
    }
    tmp_new_clusters.clear();
  } // closes loop over elements of remain used to determine connected components of remain
  return remain_clusters;
}



void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& left_new_clusters, const arma::mat &A_block){
  std::vector<std::vector<int> > tmp_new_clusters; // temporarily used to build the connected components.
  tmp_new_clusters.push_back(std::vector<int>(1, element));
  arma::mat A_tmp; 
  // loop over the existing connected components
  for(int cc = 0; cc < left_new_clusters.size(); cc++){
    A_tmp = Submatrix(A_block, tmp_new_clusters[0].size(), left_new_clusters[cc].size(), tmp_new_clusters[0], left_new_clusters[cc]);
    if(any(vectorise(A_tmp) == 1.0)){ // something in left_new_clusters[cc] is adjacent to element so we should combine these clusters
      for(int ix = 0; ix < left_new_clusters[cc].size(); ix++){
        tmp_new_clusters[0].push_back(left_new_clusters[cc][ix]);
      }
    } else{ // left_new_clusters[cc] remains its own distinct component
      tmp_new_clusters.push_back(left_new_clusters[cc]);
    }
    left_new_clusters[cc].clear();
  } // closes loop over elements of left_new_clusters
  left_new_clusters.clear();
  // update left_new_clusters: copy tmp_new_clusters
  for(int cc = 0; cc < tmp_new_clusters.size(); cc++){
    left_new_clusters.push_back(tmp_new_clusters[cc]);
    tmp_new_clusters[cc].clear();
  }
  tmp_new_clusters.clear();
  
  return ;
}

// U is an n x d matrix, where rows represent observations
// means is a d x num_splits matrix, where columns represent cluster centroids
// final_status == false means that all of our runs of kmeans failed. final_status = true means we can proceed

// outside of kmeans_repeat, we can actually return the clusters.
// really we just need it to return means.

void kmeans_repeat(arma::mat &U, arma::mat &means, bool &final_status,int num_splits,int reps){
  int n = U.n_rows; // each row of U is an observation. note that arma::kmeans wants the observations stored as columns.
  int d = U.n_cols; // dimension of the data
  arma::mat tmp_means(d, num_splits); // arma::kmeans returns the centroids stored as column vectors
  bool status = true;
  final_status = false; // if it returns false,
  arma::mat cluster_dist = arma::zeros<arma::vec>(n); // stores distance of each observation to each cluster centroid
  arma::uvec cluster_dist_indices(n); // used to sort the distance from each point to each cluster centroid
  
  double score = 0.0;
  double min_score = 0.0;
  
  for(int r = 0; r < reps; r++){
    score = 0.0;
    status = arma::kmeans(tmp_means, U.t(), num_splits, arma::random_subset, 10, false);
    if(status == false){
      Rcpp::Rcout << "kmeans failed";
    } else{
      final_status = true;
      for(int i = 0; i < n; i++){
        cluster_dist.zeros();
        for(int k = 0; k < num_splits; k++){
          cluster_dist(k) = arma::norm(U.row(i).t() - means.col(k));
        }
        cluster_dist_indices = arma::sort_index(cluster_dist, "ascend");
        score += cluster_dist(cluster_dist_indices(0)) * cluster_dist(cluster_dist_indices(0)); // adds sq distance from U(i,) to its centroid to running sum
      }
      if(r == 0){
        min_score = score;
        means = tmp_means;
      }
      else if(score < min_score){
        min_score = score;
        means = tmp_means;
      }
    }
  }
}

void kmeans_plus_plus(arma::mat &U, arma::mat &means, bool &final_status, int num_splits, int reps){
  int n = U.n_rows;
  int d = U.n_cols;
 // dimension of the data
  arma::mat tmp_means(d, num_splits); // arma::kmeans returns the centroids stored as column vectors
  bool status = true;
  final_status = false; // if it returns false,
  arma::mat cluster_dist = arma::zeros<arma::vec>(n); // stores distance of each observation to each cluster centroid
  arma::uvec cluster_dist_indices(n); // used to sort the distance from each point to each cluster centroid
  
  arma::mat dist2_mat = arma::mat(n,n);
  for(int i = 0; i < n; i++){
    for(int j = i; j < n; j++){
      dist2_mat(i,j) = arma::norm(U.row(i) - U.row(j)) * arma::norm(U.row(i) - U.row(j));
    }
  }
  arma::vec cum_prob = arma::zeros<arma::vec>(n);
  arma::vec tmp_u = arma::randu<arma::vec>(1);
  std::vector<int> already_included(n,0);
  int starting_index = 0;
  double score = 0.0;
  double min_score = 0.0;
  for(int r = 0; r < reps; r++){
    // draw a center
    tmp_u = arma::randu<arma::vec>(1);
    starting_index = floor(tmp_u(0)*n);
    already_included[starting_index] = 1;
    // starting point is indexed by floor(n*tmp_u(0))
    means.col(0) = U.row(starting_index).t();
    // now get the cumulative probabilities
    cum_prob = arma::cumsum(dist2_mat.col(starting_index))/arma::accu(dist2_mat.col(starting_index));
    for(int ns = 1; ns < num_splits; ns++){
      tmp_u = arma::randu<arma::vec>(1);
      // find first index n where tmp_u > cum_prob[n]. then use n-1
    }
  } // closes loop over r
}








