#include <thread>
#include "matlab_util.h"
#include "sdp_branch_and_bound.h"
#include "sdp_solver_util.h"
#include "JobQueue.h"
#include "Kmeans.h"
#include "util.h"
#include "config_params.h"
#include "Node.h"
#include "ThreadPool.h"

// root
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory,
                    arma::mat &C, arma::sp_mat &A, arma::vec &b, int k, int original_n, double original_trace, double global_ub) {

    // convert data
    matlab::data::TypedArray<double> C_matlab = arma_to_matlab_matrix(factory, C);
    matlab::data::CellArray A_matlab = arma_to_matlab_cell(factory, A);
    matlab::data::TypedArray<double> b_matlab = arma_to_matlab_vector(factory, b);

    // Pass vector containing args in vector
    std::vector<matlab::data::Array> args({
        factory.createScalar<int>(sdp_solver_session_threads_root),
        C_matlab, A_matlab, b_matlab,
        factory.createScalar<double>(k),
        factory.createScalar<double>(original_n),
        factory.createScalar<double>(original_trace),
        factory.createScalar<double>(sdp_solver_tol),
        factory.createScalar<int>(sdp_solver_verbose),
        factory.createScalar<int>(sdp_solver_max_cp_iter_root),
        factory.createScalar<double>(sdp_solver_cp_tol),
        factory.createScalar<double>(global_ub),
        factory.createScalar<double>(branch_and_bound_tol),
        factory.createScalar<double>(sdp_solver_eps_ineq),
        factory.createScalar<double>(sdp_solver_eps_active),
        factory.createScalar<int>(sdp_solver_max_ineq),
        factory.createScalar<int>(sdp_solver_max_pair_ineq),
        factory.createScalar<double>(sdp_solver_pair_perc),
        factory.createScalar<int>(sdp_solver_max_triangle_ineq),
        factory.createScalar<double>(sdp_solver_triangle_perc)});

    // Call MATLAB function and return result
    const size_t n_return = 11;
    matlabPtr->eval(u"clear");
    std::vector<matlab::data::Array> result = matlabPtr->feval(u"solve_cluster_cp", n_return, args);

    matlab::data::Array bound =  result[0];
    matlab::data::TypedArray<double> X_matlab = result[1];
    matlab::data::Array flag =  result[2];
    matlab::data::Array ineq = result[3];
    matlab::data::Array iter = result[4];
    matlab::data::Array iter_flag = result[5];
    matlab::data::Array pair = result[6];
    matlab::data::Array triangle = result[7];
    matlab::data::Array clique = result[8];
    matlab::data::CellArray B_matlab = result[9];
    matlab::data::TypedArray<double> l_matlab = result[10];

    double lower_bound = (double) bound[0];
    arma::mat X = matlab_to_arma_matrix(X_matlab);
    int info_flag = (int) flag[0];
    int cp_iter = (int) iter[0];
    int cp_flag = (int) iter_flag[0];
    double n_pair = (double) pair[0];
    double n_triangle = (double) triangle[0];
    double n_clique = (double) clique[0];
    int n_ineq = (int) ineq[0];
    std::vector<arma::sp_mat> B_vector = matlab_to_arma_sp_mat_vector(B_matlab);
    arma::vec l_vec = matlab_to_arma_vector(l_matlab);
    return SDPResult{info_flag, X, lower_bound, n_ineq, cp_iter, cp_flag, n_pair, n_triangle, n_clique, B_vector, l_vec};
}

// lower bound must link and cannot link
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory,
        arma::mat &C, arma::sp_mat &A, arma::vec &b, int k, int original_n, double original_trace, double global_ub,
        std::vector<arma::sp_mat> &parent_B_vector, arma::vec &parent_l_vec, int parent_n, int i, int j) {

    // convert data
    matlab::data::TypedArray<double> C_matlab = arma_to_matlab_matrix(factory, C);
    matlab::data::CellArray A_matlab = arma_to_matlab_cell(factory, A);
    matlab::data::TypedArray<double> b_matlab = arma_to_matlab_vector(factory, b);
    matlab::data::CellArray parent_Bcell = arma_to_matlab_cell(factory, parent_B_vector);
    matlab::data::TypedArray<double> parent_l = arma_to_matlab_vector(factory, parent_l_vec);

    // Pass vector containing args in vector
    std::vector<matlab::data::Array> args({
        factory.createScalar<int>(sdp_solver_session_threads),
        C_matlab, A_matlab, b_matlab,
        factory.createScalar<double>(k),
        factory.createScalar<double>(original_n),
        factory.createScalar<double>(original_trace),
        factory.createScalar<double>(sdp_solver_tol),
        factory.createScalar<int>(sdp_solver_verbose),
        factory.createScalar<int>(sdp_solver_max_cp_iter),
        factory.createScalar<double>(sdp_solver_cp_tol),
        factory.createScalar<double>(global_ub),
        factory.createScalar<double>(branch_and_bound_tol),
        factory.createScalar<double>(sdp_solver_eps_ineq),
        factory.createScalar<double>(sdp_solver_eps_active),
        factory.createScalar<int>(sdp_solver_max_ineq),
        factory.createScalar<int>(sdp_solver_max_pair_ineq),
        factory.createScalar<double>(sdp_solver_pair_perc),
        factory.createScalar<int>(sdp_solver_max_triangle_ineq),
        factory.createScalar<double>(sdp_solver_triangle_perc),
        parent_Bcell, parent_l,
        factory.createScalar<double>(parent_n),
        factory.createScalar<double>(sdp_solver_inherit_perc),
        factory.createScalar<int>(i),
        factory.createScalar<int>(j)});

    // Call MATLAB function and return result
    const size_t n_return = 11;
    matlabPtr->eval(u"clear");
    std::vector<matlab::data::Array> result = matlabPtr->feval(u"solve_cluster_cp_inherit", n_return, args);

    matlab::data::Array bound =  result[0];
    matlab::data::TypedArray<double> X_matlab = result[1];
    matlab::data::Array flag =  result[2];
    matlab::data::Array ineq = result[3];
    matlab::data::Array iter = result[4];
    matlab::data::Array iter_flag = result[5];
    matlab::data::Array pair = result[6];
    matlab::data::Array triangle = result[7];
    matlab::data::Array clique = result[8];
    matlab::data::CellArray B_matlab = result[9];
    matlab::data::TypedArray<double> l_matlab = result[10];

    double lower_bound = (double) bound[0];
    arma::mat X = matlab_to_arma_matrix(X_matlab);

    int info_flag = (int) flag[0];
    int cp_iter = (int) iter[0];
    int cp_flag = (int) iter_flag[0];
    double n_pair = (double) pair[0];
    double n_triangle = (double) triangle[0];
    double n_clique = (double) clique[0];
    int n_ineq = (int) ineq[0];
    std::vector<arma::sp_mat> B_vector = matlab_to_arma_sp_mat_vector(B_matlab);
    arma::vec l_vec = matlab_to_arma_vector(l_matlab);
    return SDPResult{info_flag, X, lower_bound, n_ineq, cp_iter, cp_flag, n_pair, n_triangle, n_clique, B_vector, l_vec};
}


// our sdp-based heuristic
double cluster_recovery(arma::mat &Ws, arma::mat &Ws_shr, arma::mat &Xopt, int k,
                        std::vector<std::pair<int, int>> &ml, std::vector<std::pair<int, int>> &cl,
                        arma::sp_mat &assignment_X, arma::mat &out_centroids) {

    int n = Ws_shr.n_rows;

    // spectral decomposition of Xopt
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, Xopt);


    arma::mat U = eigvec.cols(n - k, n - 1);
    arma::mat D = arma::diagmat(eigval);
    D = D(arma::span(n - k, n - 1), arma::span(n - k, n - 1));

    arma::mat new_Xopt = U * D * U.t();
    arma::mat C = new_Xopt * Ws_shr;

    // cluster the rows of C
    std::vector<std::pair<int, int>> ml_init = {};
    std::vector<std::pair<int, int>> cl_init = {};
    Kmeans kmeans_init(C, k, ml_init, cl_init, kmeans_verbose);
    kmeans_init.start(kmeans_max_iter, kmeans_n_start);
    arma::mat init_centroid = kmeans_init.getCentroids();

    // now perform constrained k-means with the smart initialization
    Kmeans kmeans(Ws, k, ml, cl, kmeans_verbose);
    bool flag_partition = kmeans.start(kmeans_max_iter, init_centroid);
    if (!flag_partition)
        return std::numeric_limits<double>::infinity();
    assignment_X = kmeans.getAssignments();
    out_centroids = kmeans.getCentroids();
    return kmeans.getLoss();
}


std::pair<JobData *, JobData *> create_cl_ml_jobs(double node_gap, SDPNode *node, arma::mat &X,
                                                  NodeData *parent, SharedData *shared_data, InputData *input_data) {

    if (node_gap <= branch_and_bound_tol) {
        // std::cout << "PRUNING " << node->id << "\n";
        delete(node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::make_pair(nullptr, nullptr);
    }

    std::pair<int, int> var = find_branch_norm(X);

    int i = var.first;
    int j = var.second;

    if (i == -1 && j == -1) {

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        log_file << "PRUNING BY OPTIMALITY " << node->id << "\n";
        if (node->lb - shared_data->global_ub <= -branch_and_bound_tol) {
            // update global upper bound, run the heuristic instead of setting global_ub = node->lb
            double new_ub = cluster_recovery(input_data->Ws, node->Ws, X, input_data->k,
                                             node->global_ml_pairs, node->global_cl_pairs,
                                             node->assignment_X, node->centroids);
            shared_data->global_ub = std::min(shared_data->global_ub, new_ub);
        }
        delete (node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::make_pair(nullptr, nullptr);

        // mutex is automatically released when lock goes out of scope

    }

    /*
    auto *copy_node_cl = new SDPNode(*node);
    auto *copy_node_ml = new SDPNode(*node);
    delete (node);
    */

    auto *cl_data = new NodeData();
    cl_data->node = new SDPNode(*node);
    cl_data->i = i;
    cl_data->j = j;

    auto *ml_data = new NodeData();
    ml_data->node = new SDPNode(*node);
    ml_data->i = i;
    ml_data->j = j;

    auto *cl_job_data = new JobData();
    cl_job_data->type = CANNOT_LINK;
    cl_job_data->node_data = cl_data;

    auto *ml_job_data = new JobData();
    ml_job_data->type = MUST_LINK;
    ml_job_data->node_data = ml_data;

    if (parent != nullptr) {
        delete (parent->node);
        delete (parent);
    }

    delete (node);

    return std::make_pair(cl_job_data, ml_job_data);

}


std::pair<JobData *, JobData *> build_cl_problem(MatlabStruct *matlab_struct, NodeData *node_data, InputData *input_data, SharedData  *shared_data) {

    // generate cannot link child
    auto cl_node = new SDPNode();
    cl_node->ml_map = node_data->node->ml_map;
    cl_node->local_cl_pairs = node_data->node->local_cl_pairs;
    cl_node->local_cl_pairs.emplace_back(node_data->i, node_data->j);
    cl_node->global_ml_pairs = build_global_must_link_pairs(cl_node->ml_map);
    cl_node->global_cl_pairs = build_global_cannot_link_pairs(cl_node->ml_map, cl_node->local_cl_pairs);
    // inherit Ws, A and add the cannot link constraint
    cl_node->Ws = node_data->node->Ws;
    cl_node->A = build_A_cannot_link(node_data->node->A, node_data->i, node_data->j);
    cl_node->b = build_b_cannot_link(node_data->node->b);

    auto start_time = std::chrono::high_resolution_clock::now();

    // cl_node->l = node_data->node->l + 1;
    arma::mat C = cl_node->Ws * cl_node->Ws.t();
    SDPResult sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
                                     C, cl_node->A, cl_node->b, input_data->k, input_data->Ws.n_rows, input_data->C_trace,
                                     shared_data->global_ub, node_data->node->B_vector, node_data->node->l_vec,
                                     node_data->node->Ws.n_rows, node_data->i, node_data->j);

    int flag = sdp_result.flag;
    double n_pair = sdp_result.n_pair;
    double n_triangle = sdp_result.n_triangle;
    double n_clique = sdp_result.n_clique;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    cl_node->lb = std::max(sdp_result.lb + input_data->C_trace, node_data->node->lb);
    cl_node->B_vector = sdp_result.B_vector;
    cl_node->l_vec = sdp_result.l_vec;

    if (kmeans_sdp_based)
        cl_node->ub = cluster_recovery(input_data->Ws, cl_node->Ws, sdp_result.X, input_data->k,
                                       cl_node->global_ml_pairs, cl_node->global_cl_pairs,
                                       cl_node->assignment_X, cl_node->centroids);
    else
        cl_node->ub = heuristic_solve(input_data->Ws, input_data->k, kmeans_verbose, kmeans_n_start, kmeans_max_iter,
                                      cl_node->global_ml_pairs, cl_node->global_cl_pairs,
                                      cl_node->assignment_X, cl_node->centroids);


    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        bool ub_updated = false;
        if (cl_node->ub - shared_data->global_ub <= -branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = cl_node->ub;
            shared_data->global_X = cl_node->assignment_X;
            ub_updated = true;
        }

        cl_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        shared_data->sum_ineq += n_pair + n_triangle + n_clique;
        shared_data->sum_cp_iter += cp_iter;

        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - cl_node->lb) / shared_data->global_ub;

        double gap = node_gap;
        Node *min_lb_node = shared_data->queue->getMinLb();
        if (min_lb_node != nullptr)
            gap = (shared_data->global_ub - min_lb_node->lb) / shared_data->global_ub;

        shared_data->gap = gap;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        double time = duration.count();

        //std::cout << '\r';
        print_log_sdp(log_file, cl_node->Ws.n_rows, node_data->node->id, cl_node->id, node_data->node->lb, cl_node->lb,
                      flag, time, cp_iter, cp_flag, n_ineq, n_pair, n_triangle, n_clique, shared_data->global_ub,
                      node_data->i, node_data->j, node_gap, shared_data->gap, open, ub_updated);

    }

    // mutex is automatically released when lock goes out of scope

    return create_cl_ml_jobs(node_gap, cl_node, sdp_result.X, node_data, shared_data, input_data);

}


std::pair<JobData *, JobData *> build_ml_problem(MatlabStruct *matlab_struct, NodeData *node_data, InputData *input_data, SharedData *shared_data) {

    // generate must link child
    SDPNode *ml_node;
    ml_node = new SDPNode();
    ml_node->Ws = build_Ws_must_link(node_data->node->Ws, node_data->i, node_data->j);
    ml_node->ml_map = build_must_link_map(node_data->node->ml_map, node_data->i, node_data->j);
    ml_node->local_cl_pairs = node_data->node->local_cl_pairs;
    std::set<int> dup_indices = update_cannot_link(ml_node->local_cl_pairs, node_data->i, node_data->j);
    ml_node->global_cl_pairs = build_global_cannot_link_pairs(ml_node->ml_map, ml_node->local_cl_pairs);
    ml_node->global_ml_pairs = build_global_must_link_pairs(ml_node->ml_map);
    int n = node_data->node->Ws.n_rows;
    if (n - 1 == input_data->k) {
        // build directly the assignment matrix
        // ml_node->l = node_data->node->l + 1;
        ml_node->ub = build_X_from_ml(input_data->Ws, ml_node->ml_map, ml_node->assignment_X);
        ml_node->lb = std::numeric_limits<double>::infinity();

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        ml_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;

        int open = shared_data->queue->getSize();

        bool ub_updated = false;
        if (ml_node->ub - shared_data->global_ub <= -branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = ml_node->ub;
            shared_data->global_X = ml_node->assignment_X;
            ub_updated = true;
        }

        double node_gap = (shared_data->global_ub - ml_node->lb) / shared_data->global_ub;

        double gap = node_gap;
        Node *min_lb_node = shared_data->queue->getMinLb();
        if (min_lb_node != nullptr)
            gap = (shared_data->global_ub - min_lb_node->lb) / shared_data->global_ub;

        shared_data->gap = gap;

        print_log_sdp(log_file, ml_node->Ws.n_rows, node_data->node->id, ml_node->id, node_data->node->lb,
                      ml_node->lb,0, 0, 0, 0, 0, 0, 0, 0,
                      shared_data->global_ub, node_data->i, node_data->j, node_gap, shared_data->gap, open, ub_updated);

        delete (ml_node);
        delete (node_data->node);
        return std::make_pair(nullptr, nullptr);

    }

    arma::sp_mat TTt = build_TTt(ml_node->ml_map);
    ml_node->A = build_A_must_link(n, TTt);
    ml_node->b = build_b_must_link(n, input_data->k);
    for (auto &elem : ml_node->local_cl_pairs) {
        ml_node->A = build_A_cannot_link(ml_node->A, elem.first, elem.second);
        ml_node->b = build_b_cannot_link(ml_node->b);
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // ml_node->l = node_data->node->l + 1;
    arma::mat C = ml_node->Ws * ml_node->Ws.t();

    SDPResult sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
                                     C, ml_node->A, ml_node->b, input_data->k, input_data->Ws.n_rows, input_data->C_trace,
                                     shared_data->global_ub,node_data->node->B_vector, node_data->node->l_vec,
                                     node_data->node->Ws.n_rows, node_data->i, node_data->j);

    int flag = sdp_result.flag;
    double n_pair = sdp_result.n_pair;
    double n_triangle = sdp_result.n_triangle;
    double n_clique = sdp_result.n_clique;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    ml_node->lb = std::max(sdp_result.lb + input_data->C_trace, node_data->node->lb);
    ml_node->B_vector = sdp_result.B_vector;
    ml_node->l_vec = sdp_result.l_vec;

    if (kmeans_sdp_based)
        ml_node->ub = cluster_recovery(input_data->Ws, ml_node->Ws, sdp_result.X, input_data->k,
                                       ml_node->global_ml_pairs, ml_node->global_cl_pairs,
                                       ml_node->assignment_X, ml_node->centroids);
    else
        ml_node->ub = heuristic_solve(input_data->Ws, input_data->k, kmeans_verbose, kmeans_n_start, kmeans_max_iter,
                                      ml_node->global_ml_pairs, ml_node->global_cl_pairs,
                                      ml_node->assignment_X, ml_node->centroids);


    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        bool ub_updated = false;
        if (ml_node->ub - shared_data->global_ub <= -branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = ml_node->ub;
            shared_data->global_X = ml_node->assignment_X;
            ub_updated = true;
        }

        ml_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        shared_data->sum_ineq += n_pair + n_triangle + n_clique;
        shared_data->sum_cp_iter += cp_iter;

        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - ml_node->lb) / shared_data->global_ub;

        double gap = node_gap;
        Node *min_lb_node = shared_data->queue->getMinLb();
        if (min_lb_node != nullptr)
            gap = (shared_data->global_ub - min_lb_node->lb) / shared_data->global_ub;

        shared_data->gap = gap;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        double time = duration.count();

        //std::cout << '\r';
        print_log_sdp(log_file, ml_node->Ws.n_rows, node_data->node->id, ml_node->id, node_data->node->lb, ml_node->lb,
                      flag, time, cp_iter, cp_flag, n_ineq, n_pair, n_triangle, n_clique, shared_data->global_ub,
                      node_data->i, node_data->j, node_gap, shared_data->gap, open, ub_updated);

    }

    // mutex is automatically released when lock goes out of scope

    return create_cl_ml_jobs(node_gap, ml_node, sdp_result.X, node_data, shared_data, input_data);


}

std::pair<JobData *, JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data) {

    // number of data points
    int n = input_data->Ws.n_rows;
    // init root
    SDPNode *root;
    root = new SDPNode();
    root->id = shared_data->n_nodes;
    root->Ws = input_data->Ws;
    root->A = build_A(n + 1, n);
    root->b = build_b(n + 1, input_data->k);
    for (int i = 0; i < n; i++)
        root->ml_map.insert(std::pair<int, std::set<int>> (i, {i}));
    root->global_ml_pairs = {};
    root->global_cl_pairs = {};
    root->local_cl_pairs = {};
    //root->l = 0;
    arma::mat C = root->Ws * root->Ws.t();
    double C_trace = arma::trace(C);

    auto start_time = std::chrono::high_resolution_clock::now();


    root->ub = heuristic_solve(input_data->Ws, input_data->k, kmeans_verbose, kmeans_n_start, kmeans_max_iter,
                               root->global_ml_pairs, root->global_cl_pairs,
                               root->assignment_X, root->centroids);
                               

    shared_data->global_ub = root->ub;
    shared_data->global_X = root->assignment_X;


    SDPResult sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
                                     C, root->A, root->b, input_data->k, input_data->Ws.n_rows,
                                     input_data->C_trace, shared_data->global_ub);
    int flag = sdp_result.flag;
    root->lb = sdp_result.lb + C_trace;
    double n_pair = sdp_result.n_pair;
    double n_triangle = sdp_result.n_triangle;
    double n_clique = sdp_result.n_clique;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    shared_data->n_nodes++;
    shared_data->sum_cp_iter += cp_iter;
    shared_data->sum_ineq += n_pair + n_triangle + n_clique;
    root->B_vector = sdp_result.B_vector;
    root->l_vec  = sdp_result.l_vec;

    if (kmeans_sdp_based)
        root->ub = cluster_recovery(input_data->Ws, root->Ws, sdp_result.X, input_data->k,
                                    root->global_ml_pairs, root->global_cl_pairs,
                                    root->assignment_X, root->centroids);
    else
        root->ub = heuristic_solve(input_data->Ws, input_data->k, kmeans_verbose, kmeans_n_start, kmeans_max_iter,
                                   root->global_ml_pairs, root->global_cl_pairs,
                                   root->assignment_X, root->centroids);


    if (root->ub - shared_data->global_ub <= -branch_and_bound_tol) {
        // update global upper bound
        shared_data->global_ub = root->ub;
        shared_data->global_X = root->assignment_X;
    }

    int open = shared_data->queue->getSize();

    double node_gap = (shared_data->global_ub - root->lb) / shared_data->global_ub;
    shared_data->gap = node_gap;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    double time = duration.count();

    // std::cout << '\r';
    print_log_sdp(log_file, n, -1, root->id, -std::numeric_limits<double>::infinity(), root->lb,
                  flag, time, cp_iter, cp_flag, n_ineq, n_pair, n_triangle, n_clique, shared_data->global_ub, -1, -1,
                  node_gap, node_gap, open, true);

    return create_cl_ml_jobs(node_gap, root, sdp_result.X, nullptr, shared_data, input_data);

}

bool is_thread_pool_working(std::vector<bool> &thread_state) {
    int count = 0;
    for (auto && i : thread_state) {
        if (i)
            count++;
    }
    if (count == 0)
        return false;
    return true;
}


arma::sp_mat sdp_branch_and_bound(int k, arma::mat &Ws) {

    int n_thread = branch_and_bound_parallel;

    JobAbstractQueue *queue = nullptr;
    switch (branch_and_bound_visiting_strategy) {
        case DEPTH_FIRST:
            queue = new JobStack();
            break;
        case BEST_FIRST:
            queue = new JobPriorityQueue();
            break;
        case BREADTH_FIRST:
            queue = new JobQueue();
            break;
        default:
            queue = nullptr;
    }

    auto *shared_data = new SharedData();
    shared_data->global_ub = std::numeric_limits<double>::infinity();
    shared_data->n_nodes = 0;
    shared_data->sum_ineq = 0.0;
    shared_data->sum_cp_iter = 0.0;
    shared_data->queue = queue;

    shared_data->threadStates.reserve(n_thread);
    for (int i = 0; i < n_thread; i++) {
        shared_data->threadStates.push_back(false);
    }
    
    arma::mat C = Ws * Ws.t();
    double C_trace = arma::trace(C);

    auto *input_data = new InputData();
    input_data->Ws = Ws;
    input_data->C_trace = C_trace;
    input_data->k = k;

    ThreadPool pool(shared_data, input_data, n_thread);
    
    print_header_sdp(log_file);

    auto start_all = std::chrono::high_resolution_clock::now();
    
    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(sdp_solver_folder);
    
    std::pair<JobData *, JobData *> jobs = build_root_problem(matlab_struct, input_data, shared_data);

    delete (matlab_struct);
    
    double root_gap = shared_data->gap;

    JobData *cl_job = jobs.first;
    JobData *ml_job = jobs.second;
    if (cl_job != nullptr && ml_job != nullptr) {
        pool.addJob(cl_job);
        pool.addJob(ml_job);
    }

    while (true) {

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (is_thread_pool_working(shared_data->threadStates) && shared_data->n_nodes < branch_and_bound_max_nodes) {
                shared_data->mainConditionVariable.wait(l);
            }

            if (shared_data->queue->empty() || shared_data->n_nodes >= branch_and_bound_max_nodes)
                break;
        }

    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto duration_all = std::chrono::duration_cast<std::chrono::seconds>(end_all - start_all);

    pool.quitPool();

    if (queue->empty())
        shared_data->gap = 0.0;

    log_file << "\n";
    log_file << "WALL_TIME: " << duration_all.count() << " sec\n";
    log_file << "N_NODES: " << shared_data->n_nodes << "\n";
    log_file << "AVG_INEQ: " << (double) shared_data->sum_ineq / shared_data->n_nodes << "\n";
    log_file << "AVG_CP_ITER: " << (double) shared_data->sum_cp_iter / shared_data->n_nodes << "\n";
    log_file << "ROOT_GAP: " << std::max(0.0, root_gap) << "\n";
    log_file << "GAP: " << std::max(0.0, shared_data->gap) << "\n";
    log_file << "BEST: " << shared_data->global_ub << "\n\n";

    arma::sp_mat result = shared_data->global_X;

    log_file << result;

    // free memory

    delete (input_data);
    delete (queue);
    delete (shared_data);

    return result;

}
