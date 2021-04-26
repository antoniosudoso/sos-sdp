#ifndef CLUSTERING_UTIL_H
#define CLUSTERING_UTIL_H

#include <map>
#include <set>
#include <vector>
#include <armadillo>
#include "Kmeans.h"

std::map<int, std::set<int>> build_must_link_map(std::map<int, std::set<int>> &old_ml_map, int i, int j);
std::vector<std::pair<int, int>> build_global_cannot_link_pairs(std::map<int, std::set<int>> &ml_map, std::vector<std::pair<int, int>> &local_cl_pairs);
std::set<int> update_cannot_link(std::vector<std::pair<int, int>> &local_cl_pairs, int i, int j);
std::vector<std::pair<int, int>> build_global_must_link_pairs(std::map<int, std::set<int>> &ml_map);
double build_X_from_ml(arma::mat &Ws, std::map<int, std::set<int>> &ml_map, arma::sp_mat &result_X);
arma::mat build_Ws_must_link(arma::mat &old_Ws, int i, int j);
arma::sp_mat build_TTt(std::map<int, std::set<int>> &ml_map);
std::pair<int, int> find_branch_abs(arma::mat &Z);
std::pair<int, int> find_branch(arma::mat &Z);
std::pair<int, int> find_branch_norm(arma::mat &Z);

double heuristic_solve(arma::mat &data, int k, bool verbose, int n_start, int max_iter,
                       std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl,
                       arma::sp_mat &assignment_X, arma::mat &out_centroids);

void print_header_sdp(std::ostream &log_file);
void print_log_sdp(std::ostream &log_file, int n, int node_parent, int node, double lb_parent, double lb,
                   int flag, double time, int cp_iter, int cp_flag, int n_ineq, double n_pair, double n_triangle, double n_clique,
                   double gub, int i, int j, double node_gap, double gap, int open, bool update);

#endif //CLUSTERING_UTIL_H
