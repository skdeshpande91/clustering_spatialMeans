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

