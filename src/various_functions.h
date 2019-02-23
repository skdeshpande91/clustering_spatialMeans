//
//  various_functions.hpp
//  
//
//  Created by Sameer Deshpande on 10/29/18.
//

#ifndef various_functions_h
#define various_functions_h
#include <RcppArmadillo.h>

#include <stdio.h>
arma::mat Submatrix(arma::mat M, int n_rows, int n_cols, std::vector<int> row_index, std::vector<int> col_index);
void Connected_Components(const arma::mat M, const int n, int* components, int* count);
void DFSUtil(const arma::mat M, const int n, int v, bool* visited, int* components, int* count);
arma::mat Distance_matrix(arma::vec alpha_hat, int nObs);

void new_Connected_Components(const arma::mat &M, const int n, std::vector<int> &init_components, std::vector<std::vector<int> >&components);
//void new_DFSUtil(arma::mat M, int n, int v, std::vector<bool> visited, std::vector<int> components, int count);
std::vector<std::vector<int> > Alternative_Connected_Components(std::vector<int> remain, const arma::mat &A_block);
void Alternative_Connected_Components(int element, std::vector<std::vector<int> >& left_new_clusters, const arma::mat &A_block);


void kmeans_repeat(arma::mat &U, arma::mat &means, bool &final_status, int num_splits, int reps);


#endif /* various_functions_hpp */
