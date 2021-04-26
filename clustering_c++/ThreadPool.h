#ifndef CLUSTERING_THREADPOOL_H
#define CLUSTERING_THREADPOOL_H

#include "sdp_branch_and_bound.h"
#include "util.h"
#include "matlab_util.h"

class ThreadPool {

private:

    InputData *input_data;
    SharedData  *shared_data;

    // We store the threads in a vector, so we can later stop them gracefully
    std::vector<std::thread> threads;

    // This will be set to true when the thread pool is shutting down. This tells
    // the threads to stop looping and finish
    bool done;

    void doWork(int id);


public:

    ThreadPool(SharedData *shared_data, InputData *input_data, int n_thread);
    void quitPool();
    void addJob(JobData *job_data);

};


#endif //CLUSTERING_THREADPOOL_H
