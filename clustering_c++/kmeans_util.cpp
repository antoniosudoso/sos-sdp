#include "kmeans_util.h"

void print_pairs(std::vector<std::pair<int, int>> &cl_vector) {
    for (auto &elem : cl_vector) {
        std::cout << "(" << elem.first << " " << elem.second << ")" << " ";
    }
    std::cout << "\n";
}

// sort the vector elements by second element of pairs
bool sort_by_value(const std::pair<int, double> &a, const std::pair<int, double> &b) {
    return (a.second < b.second);
}

// compute the l2-norm
double squared_distance(const arma::vec &a, const arma::vec &b) {
    double norm = arma::norm(a - b, 2);
    return std::pow(norm, 2);
}


void add_both(std::map<int, std::set<int>> &graph, int i, int j) {
    graph[i].insert(j);
    graph[j].insert(i);
}

void dfs(int i, std::map<int, std::set<int>> &graph, std::vector<bool> &visited, std::vector<int> &component) {
    visited[i] = true;
    for (auto &j : graph[i]) {
        if (!visited[j]) {
            dfs(j, graph, visited, component);
        }
    }
    component.push_back(i);
}


LinkConstraint transitive_closure(std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl, int n) {

    std::map<int, std::set<int>> ml_graph;
    std::map<int, std::set<int>> cl_graph;

    for (int i = 0; i < n; i++) {
        ml_graph.insert(std::pair<int, std::set<int>> (i, {}));
        cl_graph.insert(std::pair<int, std::set<int>> (i, {}));
    }

    for (auto &pair_ml : ml) {
        add_both(ml_graph, pair_ml.first, pair_ml.second);
    }

    std::vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            std::vector<int> component;
            dfs(i, ml_graph, visited, component);
            for (auto &x1 : component) {
                for (auto &x2 : component) {
                    if (x1 != x2) {
                        ml_graph[x1].insert(x2);
                    }
                }
            }
        }
    }

    for (auto &pair_cl : cl) {
        int i = pair_cl.first;
        int j = pair_cl.second;
        add_both(cl_graph, i, j);
        for (auto &y : ml_graph[j]) {
            add_both(cl_graph, i, y);
        }
        for (auto &x : ml_graph[i]) {
            add_both(cl_graph, x, j);
            for (auto &y : ml_graph[j]) {
                add_both(cl_graph, x, y);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (auto &j : ml_graph[i]) {
            std::set<int> set_j = cl_graph[i];
            if (j != i) {
                if (set_j.count(j)) {
                    std::fprintf(stderr, "Inconsistent constraints between %d and %d", i, j);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    return LinkConstraint{ml_graph, cl_graph};

}

void display_graph(std::map<int, std::set<int>> &map) {

    for (auto &map_elem : map) {
        int key = map_elem.first;
        std::set<int> value = map_elem.second;
        if (value.empty())
            continue;
        std::printf("%d: ", key);
        std::printf("{");
        for (auto &set_elem : value) {
            std::printf(" %d ", set_elem);
        }
        std::printf("}\n");
    }
}