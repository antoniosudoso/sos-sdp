#ifndef CLUSTERING_NODE_H
#define CLUSTERING_NODE_H

#include <armadillo>
#include <map>
#include <set>
#include <vector>

class Node {

public:

    // history
    std::map<int, std::set<int>> ml_map;

    // updated cannot link constraints
    std::vector<std::pair<int, int>> local_cl_pairs;

    // must link constraints for kmeans
    std::vector<std::pair<int, int>> global_ml_pairs;
    // cannot link constraints for kmeans
    std::vector<std::pair<int, int>> global_cl_pairs;

    // local data
    arma::mat Ws;

    // kmeans
    arma::sp_mat assignment_X;
    arma::mat centroids;

    // lower bound
    double lb;

    // upper bound
    double ub;

    // depth
    // int l;

    // node id
    int id;

};


class SDPNode : public Node {


public:

    arma::sp_mat A;
    arma::vec b;

    std::vector<arma::sp_mat> B_vector;
    arma::vec l_vec;

};

typedef struct NodeData {

    SDPNode *node;
    int i;
    int j;

} NodeData;

typedef struct JobData {

    int type;
    NodeData *node_data;

} JobData;


#endif //CLUSTERING_NODE_H
