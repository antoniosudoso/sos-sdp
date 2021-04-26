#ifndef CLUSTERING_CONFIG_PARAMS_H
#define CLUSTERING_CONFIG_PARAMS_H

#define ROOT -1
#define CANNOT_LINK 0
#define MUST_LINK 1

#define DEFAULT 0
#define ABSOLUTE 1
#define NORM 2

#define BEST_FIRST 0
#define DEPTH_FIRST 1
#define BREADTH_FIRST 2

// data full path
extern const char *data_path;
extern const char *log_path;
extern std::ofstream log_file;

// branch and bound
extern double branch_and_bound_tol;
extern int branch_and_bound_parallel;
extern int branch_and_bound_max_nodes;
extern int branch_and_bound_visiting_strategy;

// sdp solver
// extern const char *sdp_solver_matlab_session;
extern int sdp_solver_session_threads_root;
extern int sdp_solver_session_threads;
extern const char *sdp_solver_folder;
extern double sdp_solver_tol;
extern int sdp_solver_verbose;
// extern int sdp_solver_type;
extern int sdp_solver_max_cp_iter_root;
extern int sdp_solver_max_cp_iter;
extern double sdp_solver_cp_tol;
extern double sdp_solver_eps_ineq;
extern double sdp_solver_eps_active;
extern int sdp_solver_max_ineq;
extern int sdp_solver_max_pair_ineq;
extern double sdp_solver_pair_perc;
extern int sdp_solver_max_triangle_ineq;
extern double sdp_solver_triangle_perc;
extern double sdp_solver_inherit_perc;

// kmeans
extern bool kmeans_sdp_based;
extern int kmeans_max_iter;
extern int kmeans_n_start;
extern bool kmeans_verbose;

#endif //CLUSTERING_CONFIG_PARAMS_H
