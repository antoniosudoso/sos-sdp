#include "JobQueue.h"


bool JobAbstractQueue::compare_node(JobData *a, JobData *b) {
    return a->node_data->node->lb < b->node_data->node->lb;
}

JobData *JobAbstractQueue::pop() {
    JobData *node = queue.front();
    queue.pop_front();
    return node;
}

bool JobAbstractQueue::empty() {
    return queue.empty();
}

void JobAbstractQueue::sort() {
    std::sort(queue.begin(), queue.end(), compare_node);
}

Node *JobAbstractQueue::getMinLb() {
    double min_lb = std::numeric_limits<double>::infinity();
    Node *min_node = nullptr;
    for (auto current : queue) {
        if (current->node_data->node->lb < min_lb) {
            min_lb = current->node_data->node->lb;
            min_node = current->node_data->node;
        }
    }
    return min_node;
}

Node *JobAbstractQueue::getMaxLb() {
    double max_lb = - std::numeric_limits<double>::infinity();
    Node *max_node = nullptr;
    for (auto current : queue) {
        if (current->node_data->node->lb > max_lb) {
            max_lb = current->node_data->node->lb;
            max_node = current->node_data->node;
        }
    }
    return max_node;
}

void JobAbstractQueue::print() {
    for (auto &elem : queue) {
        std::cout << elem->node_data->node->lb << " ";
    }
    std::cout << "\n";
}

int JobAbstractQueue::getSize() {
    return queue.size();
}

void JobQueue::push(JobData *node) {
    queue.push_back(node);
}

void JobStack::push(JobData *node) {
    queue.push_front(node);
}

void JobPriorityQueue::push(JobData *node) {
    queue.push_back(node);
    sort();
}
