#ifndef CLUSTERING_SDP_BRANCH_AND_BOUND_H
#define CLUSTERING_SDP_BRANCH_AND_BOUND_H

#include <armadillo>
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include "JobQueue.h"

typedef struct MatlabStruct {

    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;

} MatlabStruct;


typedef struct SharedData {

    // Between workers and main
    std::condition_variable mainConditionVariable;
    std::vector<bool> threadStates;

    // Queue of requests waiting to be processed
    JobAbstractQueue *queue;
    // This condition variable is used for the threads to wait until there is work to do
    std::condition_variable queueConditionVariable;
    // Mutex to protect queue
    std::mutex queueMutex;

    double global_ub;
    double gap;
    arma::sp_mat global_X;
    int n_nodes;
    double sum_ineq;
    double sum_cp_iter;

} SharedData;

typedef struct InputData {

    arma::mat Ws;
    int k;
    double C_trace;

} InputData;


arma::sp_mat sdp_branch_and_bound(int k, arma::mat &Ws);
std::pair<JobData *, JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data);
std::pair<JobData *, JobData *> build_cl_problem(MatlabStruct *matlab_struct, NodeData *job_data, InputData *input_data, SharedData *shared_data);
std::pair<JobData *, JobData *> build_ml_problem(MatlabStruct *matlab_struct, NodeData *job_data, InputData *input_data, SharedData *shared_data);
#endif //CLUSTERING_SDP_BRANCH_AND_BOUND_H
