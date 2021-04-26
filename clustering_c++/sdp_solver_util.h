#ifndef CLUSTERING_SDP_SOLVER_UTIL_H
#define CLUSTERING_SDP_SOLVER_UTIL_H

#include <armadillo>
#include <set>
#include <MatlabEngine.hpp>

typedef struct SDPResult {

    int flag;
    arma::mat X;
    double lb;
    int n_ineq;
    int cp_iter;
    int cp_flag;
    double n_pair;
    double n_triangle;
    double n_clique;
    std::vector<arma::sp_mat> B_vector;
    arma::vec l_vec;

} SDPResult;

arma::sp_mat build_A(int m, int n);
arma::vec build_b(int m, int k);

arma::sp_mat build_TTt(std::map<int, std::set<int>> &ml_map);
arma::sp_mat build_A_must_link(int m, arma::sp_mat &TTt);
arma::vec build_b_must_link(int m, int k);
arma::sp_mat build_A_cannot_link(arma::sp_mat &old_A, int i, int j);
arma::vec build_b_cannot_link(arma::vec &old_b);

#endif //CLUSTERING_SDP_SOLVER_UTIL_H
