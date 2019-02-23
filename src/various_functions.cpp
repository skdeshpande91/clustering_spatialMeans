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

// new_alpha_bar: new alphabar for the newly created clusters
// orig_alpha_bar: the alphabars for all original clusters
// this function will find nearest neighbors and also sort
void get_subcluster_neighbors(const arma::mat &A_block, const int split_k,
                              LPPartition gamma_l, const double a1 = 1.0, const double a2 = 1.0,
                              const double nu_sigma = 1, const double lambda_sigma = 1,
                              const double rho = 0.99, const double lambda = 1.0,
                              const double eta = 1.0,
                              std::vector<std::vector <int> > init_new_clusters, std::vector<std::vector<int> > new_clusters)
{
  int num_new_clusters = init_new_clusters.size();
  int orig_K = orig_alpha_bar.size() ; // how many clusters originally
  new_clusters.clear();
  new_clusters.resize(num_new_clusters);
  std::vector<int> k_star(num_new_clusters, -1); // holds the id of nearest neighbor to each subcluster
  std::vector<int> tmp_nn; // holds potential nearest neighbors
  std::vector<double> tmp_dist; // holds distance to potential nearest neighbor
  arma::vec tmp_dist_vec = arma::zeros<vec>(1); // for sorting distances to nearest neighbors
  arma::uvec tmp_dist_indices(1); // for getting the index after sorting
  
  new_clusters.clear();
  new_clusters.resize(init_new_clusters.size());
  std::vector<int> k_star(init_new_clusters.size(), -1);
  
  double tmp_alpha_bar = 0.0;
  for(int new_k = 0; new_k < num_new_clusters; new_k++){
    tmp_alpha_bar = new_alpha_bar[new_k]; // alpha-bar for the new sub-cluster
    tmp_nn.clear();
    tmp_dist.clear();
    for(int kk = 0; kk < orig_K; kk++){
      if(kk != split_k){
        A_tmp = Submatrix(A_block, )
      }
    }
  }
  
  
  // need to loop over new_clusters
  // then need to find adjacent

  
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













