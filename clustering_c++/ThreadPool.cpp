#include "ThreadPool.h"
#include <thread>
#include <condition_variable>
#include <vector>
#include "ThreadPool.h"
#include "config_params.h"


ThreadPool::ThreadPool(SharedData *shared_data, InputData *input_data, int n_thread) {

    // This returns the number of threads supported by the system.
    // auto numberOfThreads = std::thread::hardware_concurrency();
    int numberOfThreads = n_thread;

    this->shared_data = shared_data;
    this->input_data = input_data;

    done = false;

    for (int i = 0; i < numberOfThreads; ++i) {
        // The threads will execute the private member `doWork`. Note that we need
        // to pass a reference to the function (namespaced with the class name) as
        // the first argument, and the current object as second argument
        threads.emplace_back(&ThreadPool::doWork, this, i);
    }

}

// The destructor joins all the threads so the program can exit gracefully.
void ThreadPool::quitPool() {

    {
        std::lock_guard<std::mutex> l(shared_data->queueMutex);
        // So threads know it's time to shut down
        done = true;
    }

    // Wake up all the threads, so they can finish and be joined
    shared_data->queueConditionVariable.notify_all();

    for (auto& thread : threads) {
        thread.join();
    }
}

// This function will be called by the server every time there is a request
// that needs to be processed by the thread pool
void ThreadPool::addJob(JobData *job) {
    // Grab the mutex
    std::lock_guard<std::mutex> l(shared_data->queueMutex);

    // Push the request to the queue
    shared_data->queue->push(job);

    // Notify one thread that there are requests to process
    shared_data->queueConditionVariable.notify_one();
}

// Function used by the threads to grab work from the queue
void ThreadPool::doWork(int id) {

    while (true) {

        JobData *job;

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (shared_data->queue->empty() && !done) {
                // Only wake up if there are elements in the queue or the program is shutting down
                shared_data->queueConditionVariable.wait(l);
            }


            // If we are shutting down exit without trying to process more work
            if (done) break;


            shared_data->threadStates[id] = true;

            /*
            for (auto val : shared_data->threadStates){
                std::cout << val << " ";
            }
            std::cout << "\n";
            */

            job = shared_data->queue->pop();
        }

        std::pair<JobData *, JobData *> jobs;
        
        auto *matlab_struct = new MatlabStruct();
        matlab_struct->matlabPtr = start_matlab(sdp_solver_folder);

        switch (job->type) {
            case CANNOT_LINK:
                jobs = build_cl_problem(matlab_struct, job->node_data, input_data, shared_data);
                break;
            case MUST_LINK:
                jobs = build_ml_problem(matlab_struct, job->node_data, input_data, shared_data);
                break;
            default:
                jobs = std::make_pair(nullptr, nullptr);
                break;
        }

        delete (matlab_struct);
        delete (job);

        JobData *cl_job = jobs.first;
        JobData *ml_job = jobs.second;
        if (cl_job != nullptr && ml_job != nullptr) {
            this->addJob(cl_job);
            this->addJob(ml_job);
        }

        {
            std::lock_guard<std::mutex> l(shared_data->queueMutex);
            shared_data->threadStates[id] = false;
        }

        shared_data->mainConditionVariable.notify_one();

    }
}
