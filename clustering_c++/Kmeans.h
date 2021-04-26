#ifndef CLUSTERING_KMEANS_H
#define CLUSTERING_KMEANS_H

#include <armadillo>
#include "kmeans_util.h"

class Kmeans {

private:

    bool verbose;
    arma::mat data;
    arma::mat centroids;
    arma::vec assignments;
    int n, d, k;
    double loss;
    LinkConstraint constraint;

    double objectiveFunction();
    bool assignPoints();
    bool computeCentroids();
    void initCentroids();
    bool violateConstraint(int point_i, int cluster_j);

public:

    Kmeans(const arma::mat &data, int k, std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl, bool verbose);
    bool start(int max_iter, int n_start);
    bool start(int max_iter, arma::mat &init_centroids);
    arma::sp_mat getAssignments();
    double getLoss();
    arma::mat getCentroids();

};



#endif //CLUSTERING_KMEANS_H
