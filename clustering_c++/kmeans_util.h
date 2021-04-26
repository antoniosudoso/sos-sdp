#ifndef CLUSTERING_KMEANS_UTIL_H
#define CLUSTERING_KMEANS_UTIL_H

#include <armadillo>
#include <set>

typedef struct LinkConstraint {

    // must-link
    std::map<int, std::set<int>> ml_graph;
    // cannot-link
    std::map<int, std::set<int>> cl_graph;

} LinkConstraint;

bool sort_by_value(const std::pair<int, double> &a, const std::pair<int, double> &b);
double squared_distance(const arma::vec &a, const arma::vec &b);
LinkConstraint transitive_closure(std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl, int n);
void display_graph(std::map<int, std::set<int>> &map);
void print_pairs(std::vector<std::pair<int, int>> &cl_vector);

#endif //CLUSTERING_KMEANS_UTIL_H
