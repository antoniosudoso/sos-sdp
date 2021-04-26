#ifndef CLUSTERING_JOBQUEUE_H
#define CLUSTERING_JOBQUEUE_H


#include <deque>
#include "Node.h"

class JobAbstractQueue {

protected:
    std::deque<JobData *> queue;

private:
    static bool compare_node(JobData *a, JobData *b);

public:

    virtual void push(JobData *node) = 0;

    virtual ~JobAbstractQueue() = default;

    JobData *pop();

    bool empty();

    void sort();

    Node *getMinLb();

    Node *getMaxLb();

    int getSize();

    void print();
};

class JobQueue : public JobAbstractQueue {

public:

    void push(JobData *node) override;

};

class JobStack : public JobAbstractQueue {

public:

    void push(JobData *node) override;

};

class JobPriorityQueue : public JobAbstractQueue {

public:

    void push(JobData *node) override;

};


#endif //CLUSTERING_JOBQUEUE_H
