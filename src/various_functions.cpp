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

void Connected_Components(arma::mat M, int n, int* components, int* count)
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



void DFSUtil(arma::mat M, int n, int v, bool* visited, int* components, int* count)
{
  visited[v] = true;
  components[v] = *count;
  for(int i = 0; i < n; i++)
    if(M(v,i) == 1.0)
      if(visited[i] == false)
        DFSUtil(M, n, i, visited, components, count);
}


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

std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain){
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

















