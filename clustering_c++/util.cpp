#include "util.h"
#include <iomanip>


// build the must link map history
std::map<int, std::set<int>> build_must_link_map(std::map<int, std::set<int>> &old_ml_map, int i, int j) {

    std::map<int, std::set<int>> ml_map;

    int n = old_ml_map.size();

    // init new map
    for (int k = 0; k < n - 1; k++){
        ml_map.insert(std::pair<int, std::set<int>>(k, {}));
    }

    for (int k = 0; k < j; k++) {
        // copy keys from old map to new map
        for (auto &elem : old_ml_map[k]) {
            ml_map[k].insert(elem);
        }
    }

    for (auto &elem : old_ml_map[j]) {
        ml_map[i].insert(elem);
    }

    // shift
    for (int k = j + 1; k < n; k++) {
        for (auto &elem : old_ml_map[k]) {
            ml_map[k - 1].insert(elem);
        }
    }

    return ml_map;
}

// build global cannot link pairs from local cannot link and map history
std::vector<std::pair<int, int>> build_global_cannot_link_pairs(std::map<int, std::set<int>> &ml_map, std::vector<std::pair<int, int>> &local_cl_pairs) {

    std::vector<std::pair<int, int>> global_cl_pairs;

    for (auto &elem : local_cl_pairs) {

        int key_first = elem.first;
        int key_second = elem.second;
        std::set<int> set_first = ml_map[key_first];
        std::set<int> set_second = ml_map[key_second];
        for (auto &set_elem_1 : set_first) {
            for (auto &set_elem_2 : set_second) {
                global_cl_pairs.emplace_back(set_elem_1, set_elem_2);
            }
        }
    }

    return global_cl_pairs;
}

// update indices in the local cannot link
std::set<int> update_cannot_link(std::vector<std::pair<int, int>> &local_cl_pairs, int i, int j) {
    int size = local_cl_pairs.size();

    for (int k = 0; k < size; k++) {

        if (local_cl_pairs[k].first == j)
            local_cl_pairs[k].first = i;
        if (local_cl_pairs[k].first > j)
            local_cl_pairs[k].first -= 1;

        if (local_cl_pairs[k].second == j)
            local_cl_pairs[k].second = i;
        if (local_cl_pairs[k].second > j)
            local_cl_pairs[k].second -= 1;

        // check right order
        if (local_cl_pairs[k].first >= local_cl_pairs[k].second) {
            int dmy = local_cl_pairs[k].second;
            local_cl_pairs[k].second = local_cl_pairs[k].first;
            local_cl_pairs[k].first = dmy;
        }

    }

    std::set<int> dup_indices = {};
    for (int t = 0; t < size; t++) {
        for (int s = t + 1; s < size; s++) {
            std::pair<int, int> p1 = local_cl_pairs[t];
            std::pair<int, int> p2 = local_cl_pairs[s];
            if (p1.first == p2.first && p1.second == p2.second) {
                dup_indices.insert(s);
            }
        }
    }

    if (!dup_indices.empty()) {
        // delete indices (indices are sorted in the set)
        std::set<int>::reverse_iterator rit;
        for (rit = dup_indices.rbegin(); rit != dup_indices.rend(); rit++)
            local_cl_pairs.erase(local_cl_pairs.begin() + *rit);
    }

    return dup_indices;

}

// build global must link from history map
std::vector<std::pair<int, int>> build_global_must_link_pairs(std::map<int, std::set<int>> &ml_map) {

    std::vector<std::pair<int, int>> ml_pairs;

    for (auto &elem_map : ml_map) {
        std::set<int> set_i = elem_map.second;
        int size = set_i.size();
        if (size > 1) {
            bool is_first = true;
            int first;
            for (auto &elem_set : set_i) {
                if (is_first) {
                    first = elem_set;
                    is_first = false;
                    continue;
                }
                ml_pairs.emplace_back(first, elem_set);
            }
        }
    }

    return ml_pairs;

}


// build directly the assignment matrix from must-link constraints, return assignment matrix in result_X
double build_X_from_ml(arma::mat &Ws, std::map<int, std::set<int>> &ml_map, arma::sp_mat &result_X) {
    int n = Ws.n_rows;
    int d = Ws.n_cols;
    int k = ml_map.size();
    arma::sp_mat X(n, k);
    arma::mat centroids = arma::zeros(k, d);
    for (int i = 0; i < k; i++) {
        std::set<int> set = ml_map[i];
        for (auto &elem : set) {
            X(elem, i) = 1;
            centroids.row(i) += Ws.row(elem);
        }
        centroids.row(i) = centroids.row(i) / set.size();
    }
    result_X = X;
    arma::mat m = Ws - X * centroids;
    return arma::dot(m.as_col(), m.as_col());
}

// sum i and j and build new Ws
arma::mat build_Ws_must_link(arma::mat &old_Ws, int i, int j) {
    if (i >= j) {
        std::fprintf(stderr, "build_Ws_must_link(): branching node must be sorted!");
        exit(EXIT_FAILURE);
    }
    int n = old_Ws.n_rows;
    int d = old_Ws.n_cols;
    arma::mat Ws(n - 1, d);
    for (int k = 0; k < j; k++)
        Ws.row(k) = old_Ws.row(k);
    Ws.row(i) += old_Ws.row(j);
    for (int k = j; k < n - 1; k++)
        Ws.row(k) = old_Ws.row(k + 1);
    return Ws;
}


arma::sp_mat build_TTt(std::map<int, std::set<int>> &ml_map) {
    int n = ml_map.size();
    arma::sp_mat T(n, n);
    int i = 0;
    for (auto &elem : ml_map) {
        T(i, i) = elem.second.size();
        i++;
    }
    return T;
}


bool comparator_find_branch(std::pair<std::pair<int, int>, double> &a, std::pair<std::pair<int, int>, double> &b) {
    return (a.second < b.second);
}

// branching decision
std::pair<int, int> find_branch(arma::mat &Z) {

    std::vector<std::pair<std::pair<int, int>, double>> branch_data;

    int n = Z.n_rows;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            double min1 = std::min(Z(i, j), Z(i, i) - Z(i, j));
            std::pair<int, int> i_j = std::pair<int, int>(i, j);
            std::pair<std::pair<int, int>, double> data1(i_j, min1);
            double min2 = std::min(Z(i, j), Z(j, j) - Z(i, j));
            std::pair<int, int> j_i = std::pair<int, int>(j, i);
            std::pair<std::pair<int, int>, double> data2(j_i, min2);
            branch_data.push_back(data1);
            branch_data.push_back(data2);
        }
    }

    auto max_elem = std::max_element(branch_data.begin(), branch_data.end(), comparator_find_branch);
    int max_i = max_elem.base()->first.first;
    int max_j = max_elem.base()->first.second;
    double max_val = max_elem.base()->second;

    std::pair<int, int> max_pair(max_i, max_j);

    if (max_val < 1e-5) {
        // Z optimal
        return {-1, -1};
    }

    if (max_j < max_i)
        max_pair = std::pair<int, int>(max_j, max_i);

    return max_pair;
}

std::pair<int, int> find_branch_abs(arma::mat &Z) {

    std::vector<std::pair<std::pair<int, int>, double>> branch_data;

    int n = Z.n_rows;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            double min1 = std::min(Z(i, j), std::abs(Z(i, i) - Z(i, j)));
            std::pair<int, int> i_j = std::pair<int, int>(i, j);
            std::pair<std::pair<int, int>, double> data1(i_j, min1);
            double min2 = std::min(Z(i, j), std::abs(Z(j, j) - Z(i, j)));
            std::pair<int, int> j_i = std::pair<int, int>(j, i);
            std::pair<std::pair<int, int>, double> data2(j_i, min2);
            branch_data.push_back(data1);
            branch_data.push_back(data2);
        }
    }

    auto max_elem = std::max_element(branch_data.begin(), branch_data.end(), comparator_find_branch);
    int max_i = max_elem.base()->first.first;
    int max_j = max_elem.base()->first.second;
    double max_val = max_elem.base()->second;

    std::pair<int, int> max_pair(max_i, max_j);

    if (max_val < 1e-5) {
        // Z optimal
        return {-1, -1};
    }

    if (max_j < max_i)
        max_pair = std::pair<int, int>(max_j, max_i);

    return max_pair;
}

/*
std::vector<std::pair<std::pair<int, int>, double>> not_in_place_sort(std::vector<std::pair<std::pair<int, int>, double>> original) {
    std::sort(original.begin(), original.end(), comparator_find_branch);
    return original;
}
*/

// branching decision
std::pair<int, int> find_branch_norm(arma::mat &Z) {

    std::vector<std::pair<std::pair<int, int>, double>> branch_data;

    int n = Z.n_rows;
    for (int i = 0; i < n ; i++) {
        for (int j = i + 1; j < n; j++) {
            double norm = std::pow(arma::norm(Z.row(i) - Z.row(j), 2), 2);
            double first = Z(i, j);
            double min1 = std::min(first, norm);
            std::pair<int, int> i_j = std::pair<int, int>(i, j);
            std::pair<std::pair<int, int>, double> data1(i_j, min1);
            branch_data.push_back(data1);

        }
    }

    auto max_elem = std::max_element(branch_data.begin(), branch_data.end(), comparator_find_branch);
    int max_i = max_elem.base()->first.first;
    int max_j = max_elem.base()->first.second;
    double max_val = max_elem.base()->second;

    std::pair<int, int> max_pair(max_i, max_j);

    /*
    auto c = not_in_place_sort(branch_data);
    for (auto &elem : c) {
        std::cout << elem.first.first << " " << elem.first.second << " " << elem.second << "\n";
    }
    */

    //std::cout << "\nMAX_VAL: " << max_val << "\n";

    if (max_val < 1e-5) {
        // Z optimal
        return {-1, -1};
    }

    if (max_j < max_i)
        max_pair = std::pair<int, int>(max_j, max_i);

    //std::cout << "\nZ_ij = " << Z(max_pair.first, max_pair.second) << "\n";
    //std::cout << "\n" << max_pair.first << " " << max_pair.second << "\n";

    return max_pair;
}

// upper bound with constrained k-means, set result matrix as argument assignment_X and out_centroids
double heuristic_solve(arma::mat &data, int k, bool verbose, int n_start, int max_iter,
                       std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl,
                       arma::sp_mat &assignment_X, arma::mat &out_centroids) {
    // call constrained k-means
    Kmeans kmeans(data, k, ml, cl, verbose);
    bool flag_partition = kmeans.start(max_iter, n_start);
    if (!flag_partition)
        return std::numeric_limits<double>::infinity();
    assignment_X = kmeans.getAssignments();
    out_centroids = kmeans.getCentroids();
    return kmeans.getLoss();
}

void print_header_sdp(std::ostream &log_file) {

    log_file << "\n" << "|" <<
              std::setw(5) << "N" << "|" <<
              std::setw(9) << "NODE_PAR" << "|" <<
              std::setw(8) << "NODE" << "|" <<
              std::setw(12) << "LB_PAR" << "|" <<
              std::setw(12) << "LB" << "|" <<
              std::setw(6) << "FLAG" << "|" <<
              std::setw(10) << "TIME (s)" << "|" <<
              std::setw(8) << "CP_ITER" << "|" <<
              std::setw(8) << "CP_FLAG" << "|" <<
              std::setw(10) << "CP_INEQ" << "|" <<
              std::setw(9) << "PAIR" << " " <<
              std::setw(9) << "TRIANGLE" << " " <<
              std::setw(9) << "CLIQUE" << "|" <<
              std::setw(12) << "GUB" << "|" <<
              std::setw(6) << "I" << " " <<
              std::setw(6) << "J" << "|" <<
              std::setw(13) << "NODE_GAP" << "|" <<
              std::setw(13) << "GAP" << "|" <<
              std::setw(6) << "OPEN" << "|"
              << std::endl;

}

void print_log_sdp(std::ostream &log_file, int n, int node_parent, int node, double lb_parent, double lb,
                   int flag, double time, int cp_iter, int cp_flag, int n_ineq, double n_pair, double n_triangle, double n_clique,
                   double gub, int i, int j, double node_gap, double gap, int open, bool update) {

    if (!update) {

        log_file << "|" <<
                  std::setw(5) << n << "|" <<
                  std::setw(9) << node_parent << "|" <<
                  std::setw(8) << node << "|" <<
                  std::setw(12) << lb_parent << "|" <<
                  std::setw(12) << lb << "|" <<
                  std::setw(6) << flag << "|" <<
                  std::setw(10) << time << "|" <<
                  std::setw(8) << cp_iter << "|" <<
                  std::setw(8) << cp_flag << "|" <<
                  std::setw(10) << n_ineq << "|" <<
                  std::setw(9) << n_pair << " " <<
                  std::setw(9) << n_triangle << " " <<
                  std::setw(9) << n_clique << "|" <<
                  std::setw(12) << gub << "|" <<
                  std::setw(6) << i << " " <<
                  std::setw(6) << j << "|" <<
                  std::setw(13) << node_gap << "|" <<
                  std::setw(13) << gap << "|" <<
                  std::setw(6) << open << "|"
                  << std::endl;

    } else {

        log_file << "|" <<
                  std::setw(5) << n << "|" <<
                  std::setw(9) << node_parent << "|" <<
                  std::setw(8) << node << "|" <<
                  std::setw(12) << lb_parent << "|" <<
                  std::setw(12) << lb << "|" <<
                  std::setw(6) << flag << "|" <<
                  std::setw(10) << time << "|" <<
                  std::setw(8) << cp_iter << "|" <<
                  std::setw(8) << cp_flag << "|" <<
                  std::setw(10) << n_ineq << "|" <<
                  std::setw(9) << n_pair << " " <<
                  std::setw(9) << n_triangle << " " <<
                  std::setw(9) << n_clique << "|" <<
                  std::setw(11) << gub << "*|" <<
                  std::setw(6) << i << " " <<
                  std::setw(6) << j << "|" <<
                  std::setw(13) << node_gap << "|" <<
                  std::setw(13) << gap << "|" <<
                  std::setw(6) << open << "|"
                  << std::endl;

    }
}

/*
void get_original_matrix(arma::mat &X, std::vector<std::pair<int, int>> &ml) {

    std::cout << "\n";
    std::cout << "N: " << X.n_rows << "\n";
    auto rit = ml.rbegin();
    for (; rit != ml.rend(); ++rit) {
        int i = (*rit).first;
        int j = (*rit).second;
        std::cout << "\n";
        std::cout << i << " " << j << "\n";
        arma::vec row_i = X.row(i).t();
        X.insert_rows(j, row_i.t());
        arma::vec col_i = X.col(i);
        X.insert_cols(j, col_i);
        //std::cout << X << "\n";
    }
}
*/
