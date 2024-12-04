## SOS-SDP: An Exact Solver for Minimum Sum-of-Squares Clustering

<p align="center">
  <img src="https://github.com/antoniosudoso/sos-sdp/blob/main/logo.svg" width="200" height="200" />
</p>


**SOS-SDP** is an exact algorithm based on the branch-and-bound technique for solving the Minimum Sum-of-Squares Clustering (MSSC) problem described in the paper ["SOS-SDP: an Exact Solver for Minimum Sum-of-Squares Clustering"](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2022.1166). This repository contains the C++ source code, the MATLAB scripts, and the datasets used for the experiments.

> V. Piccialli, A. M. Sudoso, A. Wiegele (2022). SOS-SDP: An Exact Solver for Minimum Sum-of-Squares Clustering, **INFORMS Journal on Computing**, https://doi.org/10.1287/ijoc.2022.1166.

## Installation
**SOS-SDP** calls the semidefinite programming solver [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/) by using the [MATLAB Engine API](https://www.mathworks.com/help/matlab/calling-matlab-engine-from-cpp-programs.html) for C++. It requires the MATLAB engine library *libMatlabEngine* and the Matlab Data Array library *libMatlabDataArray*. **SOS-SDP** uses the [Armadillo](http://arma.sourceforge.net/) library to handle matrices and linear algebra operations efficiently. Before installing Armadillo, first install OpenBLAS and LAPACK along with the corresponding development files.
**SOS-SDP** implements a configurable thread pool of POSIX threads to speed up the branch-and-bound search.

Ubuntu and Debian instructions:
1) Install MATLAB (>= 2016b).
2) Install CMake, OpenBLAS, LAPACK and Armadillo:
 ```
sudo apt-get update
sudo apt-get install cmake libopenblas-dev liblapack-dev libarmadillo-dev
```
3) Open the makefile `clustering_c++/Makefile` and set the variable `matlab_path` with your MATLAB root folder.
4) Compile the code:

```
cd clustering_c++/
make
```

4) Download SDPNAL+, move the folder `clustering_matlab` containing the MATLAB source code of **SOS-SDP** in the SDPNAL+ main directory and 
set the parameter `SDP_SOLVER_FOLDER` of the configuration file accordingly. This folder and its subfolders will be automatically added to the MATLAB search path when **SOS-SDP** starts.

The code has been tested on Ubuntu Server 20.04 with MATLAB R2020b and Armadillo 10.2.

## Configuration
Various parameters used in **SOS-SDP** can be modified in the configuration file `clustering_c++/config.txt`:

- `BRANCH_AND_BOUND_TOL` - optimality tolerance of the branch-and-bound
- `BRANCH_AND_BOUND_PARALLEL` -  thread pool size: single thread (1), multi-thread (> 1)
- `BRANCH_AND_BOUND_MAX_NODES` - maximum number of nodes
- `BRANCH_AND_BOUND_VISITING_STRATEGY` - best first (0),  depth first (1), breadth first (2)
- `SDP_SOLVER_SESSION_THREADS_ROOT` - number of threads for the MATLAB session at the root
- `SDP_SOLVER_SESSION_THREADS` - number of threads for the MATLAB session for the ML and CL nodes
- `SDP_SOLVER_FOLDER` - full path of the SDPNAL+ folder
- `SDP_SOLVER_TOL` - accuracy of SDPNAL+
- `SDP_SOLVER_VERBOSE` - do not display log (0), display log (1)
- `SDP_SOLVER_MAX_CP_ITER_ROOT` - maximum number of cutting-plane iterations at the root
- `SDP_SOLVER_MAX_CP_ITER` - maximum number of cutting-plane iterations for the ML and CL nodes
- `SDP_SOLVER_CP_TOL` - cutting-plane tolerance between two consecutive cutting-plane iterations
- `SDP_SOLVER_MAX_INEQ` - maximum number of valid inequalities to add
- `SDP_SOLVER_INHERIT_PERC` - fraction of inequalities to inherit
- `SDP_SOLVER_EPS_INEQ` - tolerance for checking the violation of the inequalities
- `SDP_SOLVER_EPS_ACTIVE` - tolerance for detecting the active inequalities
- `SDP_SOLVER_MAX_PAIR_INEQ` - maximum number of pair inequalities to separate
- `SDP_SOLVER_PAIR_PERC` - fraction of the most violated pair inequalities to add
- `SDP_SOLVER_MAX_TRIANGLE_INEQ` - maximum number of triangle inequalities to separate
- `SDP_SOLVER_TRIANGLE_PERC` - fraction of the most violated triangle inequalities to add
- `KMEANS_SDP_BASED` - constrained k-means with k-means++ initialization (0), sdp-based initialization (1)
- `KMEANS_MAX_ITER` - maximum number of k-means iterations
- `KMEANS_N_START` - number of times k-means is run
- `KMEANS_VERBOSE` - do not display log (0), display log (1)

## Usage
```
cd clustering_c++/
./bb <DATASET> <K> <LOG>
```
- `DATASET` - the path of the dataset
- `K` - the number of clusters
- `LOG` - the path of the log file

The dataset file contains the data points `x_ij` and the must include an header line with the problem size `n` and the dimension `d`:

```
n d
x_11 x_12 ... x_1d
x_21 x_22 ... x_2d
...
...
x_n1 x_n2 ... x_nd
```


## Log

The log file reports the progress of the algorithm:

- `N` - size of the current node
- `NODE_PAR` - id of the parent node
- `NODE` - id of the current node
- `LB_PAR` - lower bound of the parent node
- `LB` - lower bound of the current node
- `FLAG` - termination flag of SDPNAL+
    -  `0` - SDP is solved to the required accuracy
    -  `1` - SDP is not solved successfully
    -  `-1, -2, -3` - SDP is partially solved successfully
- `TIME (s)` - computational time in seconds of the current node
- `CP_ITER` - number of cutting-plane iterations
- `CP_FLAG` - termination flag of the cutting-plane procedure
    - `-3` - current bound is worse than the previous one
    - `-2` - SDP is not solved successfully
    - `-1` - maximum number of iterations
    -  `0` - no violated inequalities
    -  `1` - maximum number of inequalities
    -  `2` - node must be pruned
    -  `3` - cutting-plane tolerance
- `CP_INEQ` - number of inequalities added in the last cutting-plane iteration
- `PAIR TRIANGLE CLIQUE` - average number of added cuts for each class of inequalities
- `GUB` - global upper bound
- `I J` - current branching decision
- `NODE_GAP` - gap at the current node
- `GAP` - overall gap 
- `OPEN` - number of open nodes

The log file prints the optimal minimum sum-of-squares (MSSC) objective and the cluster indicator matrix when the algorithm ends `(point_id, cluster_id)`.

Log file example:

```
DATA_PATH, n, d, k: /home/ubuntu/SOS-SDP/datasets/glass.txt 214 9 6
LOG_PATH: /home/ubuntu/SOS-SDP/results/glass_6.txt

BRANCH_AND_BOUND_TOL: 1e-04
BRANCH_AND_BOUND_PARALLEL: 16
BRANCH_AND_BOUND_MAX_NODES: 100
BRANCH_AND_BOUND_VISITING_STRATEGY: 0

SDP_SOLVER_SESSION_THREADS_ROOT: 16
SDP_SOLVER_SESSION_THREADS: 1
SDP_SOLVER_FOLDER: /home/ubuntu/SOS-SDP/SDPNAL+/
SDP_SOLVER_TOL: 1e-05
SDP_SOLVER_VERBOSE: 0
SDP_SOLVER_MAX_CP_ITER_ROOT: 80
SDP_SOLVER_MAX_CP_ITER: 40
SDP_SOLVER_CP_TOL: 1e-05
SDP_SOLVER_MAX_INEQ: 100000
SDP_SOLVER_INHERIT_PERC: 1
SDP_SOLVER_EPS_INEQ: 0.0001
SDP_SOLVER_EPS_ACTIVE: 1e-06
SDP_SOLVER_MAX_PAIR_INEQ: 100000
SDP_SOLVER_PAIR_PERC: 0.05
SDP_SOLVER_MAX_TRIANGLE_INEQ: 100000
SDP_SOLVER_TRIANGLE_PERC: 0.05

KMEANS_SDP_BASED: 1
KMEANS_MAX_ITER: 100
KMEANS_N_START: 50
KMEANS_VERBOSE: 0


|    N| NODE_PAR|    NODE|      LB_PAR|          LB|  FLAG|  TIME (s)| CP_ITER| CP_FLAG|   CP_INEQ|     PAIR  TRIANGLE    CLIQUE|         GUB|     I      J|     NODE_GAP|          GAP|  OPEN|
|  214|       -1|       0|        -inf|     72.9391|     0|       151|       7|       3|      3014|  639.857      4448   8.71429|    72.9709*|    -1     -1|  0.000436087|  0.000436087|     0|
|  214|        0|       1|     72.9391|     72.9644|    -1|        12|       0|       2|      3014|        0         0         0|     72.9709|   107    201|   8.8639e-05|   8.8639e-05|     0|
|  213|        0|       2|     72.9391|     72.9487|     0|        31|       1|      -3|      2485|        0       922         0|    72.9647*|   107    201|  0.000220183|  0.000220183|     0|
|  212|        2|       3|     72.9487|     72.9643|     0|         8|       0|       2|      2043|        0         0         0|     72.9647|    32     51|  5.25193e-06|  5.25193e-06|     0|
|  213|        2|       4|     72.9487|     72.9789|     0|        12|       0|       2|      2485|        0         0         0|     72.9647|    32     51|  -0.00019452|  -0.00019452|     0|

WALL_TIME: 199 sec
N_NODES: 5
AVG_INEQ: 1203.71
AVG_CP_ITER: 1.6
ROOT_GAP: 0.000436087
GAP: 0
BEST: 72.9647

```

## Related Works

> V. Piccialli, A. Russo Russo, A. M. Sudoso (2022). An exact algorithm for semi-supervised minimum sum-of-squares clustering. **Computers & Operations Research**.
- Paper: https://doi.org/10.1016/j.cor.2022.105958
- Code: https://github.com/antoniosudoso/pc-sos-sdp

> V. Piccialli, A. M. Sudoso (2023). Global optimization for cardinality-constrained minimum sum-of-squares clustering via semidefinite programming. **Mathematical Programming**.
- Paper: https://doi.org/10.1007/s10107-023-02021-8
- Code: https://github.com/antoniosudoso/cc-sos-sdp

