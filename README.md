## SOS-SDP: an Exact Solver for Minimum Sum-of-Squares Clustering

<p align="center">
  <img src="https://github.com/antoniosudoso/sos-sdp/blob/main/logo.svg" width="200" height="200" />
</p>


**SOS-SDP** is an exact algorithm based on the branch-and-bound technique for solving the Minimum Sum-of-Squares Clustering (MSSC) problem described in the paper ["SOS-SDP: an Exact Solver for Minimum Sum-of-Squares Clustering"](https://arxiv.org/abs/2104.11542). This repository contains the C++ source code, the MATLAB scripts, and the datasets used for the experiments.

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

- `BRANCH_AND_BOUND_TOL` - the optimality tolerance of the branch-and-bound
- `BRANCH_AND_BOUND_PARALLEL` -  the thread pool size: single thread (1), multi-thread (> 1)
- `BRANCH_AND_BOUND_MAX_NODES` - the maximum number of nodes
- `BRANCH_AND_BOUND_VISITING_STRATEGY` - best first (0),  depth first (1), breadth first (2)
- `SDP_SOLVER_SESSION_THREADS_ROOT` - the numer of threads for the MATLAB session at the root
- `SDP_SOLVER_SESSION_THREADS` - the numer of threads for the MATLAB session for the ML and CL nodes
- `SDP_SOLVER_FOLDER` - the full path of the SDPNAL+ folder
- `SDP_SOLVER_TOL` - the accuracy of SDPNAL+
- `SDP_SOLVER_VERBOSE` - do not display log (0), display log (1)
- `SDP_SOLVER_MAX_CP_ITER_ROOT` - the maximum number of cutting-plane iterations at the root
- `SDP_SOLVER_MAX_CP_ITER` - the maximum number of cutting-plane iterations for the ML and CL nodes
- `SDP_SOLVER_CP_TOL` - the cutting-plane tolerance between two consecutive cutting-plane iterations
- `SDP_SOLVER_MAX_INEQ` - the maximum number of valid inequalities to add
- `SDP_SOLVER_INHERIT_PERC` - the fraction of inequalities to inherit
- `SDP_SOLVER_EPS_INEQ` - the tolerance for checking the violation of the inequalities
- `SDP_SOLVER_EPS_ACTIVE` - the tolerance for detecting the active inequalities
- `SDP_SOLVER_MAX_PAIR_INEQ` - the maximum number of pair inequalities to separate
- `SDP_SOLVER_PAIR_PERC` - the fraction of the most violated pair inequalities to add
- `SDP_SOLVER_MAX_TRIANGLE_INEQ` - the maximum number of triangle inequalities to separate
- `SDP_SOLVER_TRIANGLE_PERC` - the fraction of the most violated triangle inequalities to add
- `KMEANS_SDP_BASED` - constrained k-means with k-means++ initialization (0), sdp-based initialization (1)
- `KMEANS_MAX_ITER` - the maximum number of k-means iterations
- `KMEANS_N_START` - the number of times k-means is run
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
x_11, x_12, ..., x_1d
x_21, x_22, ..., x_2d
...
...
x_n1, x_n2, ..., x_nd
```


## Log

The log file reports the progress of the algorithm:

- `N` - the size of the current node
- `NODE_PAR` - the id of the parent node
- `NODE` - the id of the current node
- `LB_PAR` - the lower bound of the parent node
- `LB` - the lower bound of the current node
- `FLAG` - the termination flag of SDPNAL+
    -  `0` - SDP is solved to the required accuracy
    -  `1` - SDP is not solved successfully
    -  `-1, -2, -3` - SDP is partially solved successfully
- `TIME (s)` - the computation time in seconds of the current node
- `CP_ITER` - the number of cutting-plane iterations
- `CP_FLAG` - the termination flag of the cutting-plane procedure
    - `-3` - current bound is worse than the previous one
    - `-2` - SDP is not solved successfully
    - `-1` - maximum number of iterations
    -  `0` - no violated inequalities
    -  `1` - maximum number of inequalities
    -  `2` - node must be pruned
    -  `3` - cutting-plane tolerance
- `CP_INEQ` - the number of inequalities added in the last cutting-plane iteration
- `PAIR TRIANGLE CLIQUE` - the average number of added cuts for each class of inequalities
- `GUB` - the global upper bound
- `I J` - the current branching decision
- `NODE_GAP` - the gap at the current node
- `GAP` - the overall gap 
- `OPEN` - the number of open nodes

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

[matrix size: 214x6; n_nonzero: 214; density: 16.7%]

    (181, 0)          1.0000
    (184, 0)          1.0000
    (190, 0)          1.0000
    (191, 0)          1.0000
    (192, 0)          1.0000
    (193, 0)          1.0000
    (194, 0)          1.0000
    (195, 0)          1.0000
    (196, 0)          1.0000
    (197, 0)          1.0000
    (198, 0)          1.0000
    (199, 0)          1.0000
    (200, 0)          1.0000
    (202, 0)          1.0000
    (203, 0)          1.0000
    (204, 0)          1.0000
    (205, 0)          1.0000
    (206, 0)          1.0000
    (207, 0)          1.0000
    (208, 0)          1.0000
    (209, 0)          1.0000
    (210, 0)          1.0000
    (211, 0)          1.0000
    (212, 0)          1.0000
    (213, 0)          1.0000
      (1, 1)          1.0000
      (2, 1)          1.0000
      (3, 1)          1.0000
      (4, 1)          1.0000
      (6, 1)          1.0000
      (7, 1)          1.0000
      (8, 1)          1.0000
      (9, 1)          1.0000
     (11, 1)          1.0000
     (14, 1)          1.0000
     (15, 1)          1.0000
     (16, 1)          1.0000
     (19, 1)          1.0000
     (22, 1)          1.0000
     (23, 1)          1.0000
     (24, 1)          1.0000
     (25, 1)          1.0000
     (26, 1)          1.0000
     (27, 1)          1.0000
     (28, 1)          1.0000
     (29, 1)          1.0000
     (31, 1)          1.0000
     (33, 1)          1.0000
     (34, 1)          1.0000
     (35, 1)          1.0000
     (37, 1)          1.0000
     (40, 1)          1.0000
     (41, 1)          1.0000
     (42, 1)          1.0000
     (45, 1)          1.0000
     (49, 1)          1.0000
     (52, 1)          1.0000
     (53, 1)          1.0000
     (54, 1)          1.0000
     (57, 1)          1.0000
     (58, 1)          1.0000
     (59, 1)          1.0000
     (60, 1)          1.0000
     (70, 1)          1.0000
     (72, 1)          1.0000
     (73, 1)          1.0000
     (74, 1)          1.0000
     (75, 1)          1.0000
     (76, 1)          1.0000
     (77, 1)          1.0000
     (79, 1)          1.0000
     (80, 1)          1.0000
     (81, 1)          1.0000
     (82, 1)          1.0000
     (83, 1)          1.0000
     (84, 1)          1.0000
     (85, 1)          1.0000
     (86, 1)          1.0000
     (87, 1)          1.0000
     (88, 1)          1.0000
     (89, 1)          1.0000
     (91, 1)          1.0000
     (93, 1)          1.0000
     (94, 1)          1.0000
     (95, 1)          1.0000
     (98, 1)          1.0000
     (99, 1)          1.0000
    (101, 1)          1.0000
    (114, 1)          1.0000
    (115, 1)          1.0000
    (116, 1)          1.0000
    (117, 1)          1.0000
    (119, 1)          1.0000
    (120, 1)          1.0000
    (122, 1)          1.0000
    (123, 1)          1.0000
    (124, 1)          1.0000
    (126, 1)          1.0000
    (132, 1)          1.0000
    (134, 1)          1.0000
    (137, 1)          1.0000
    (138, 1)          1.0000
    (139, 1)          1.0000
    (140, 1)          1.0000
    (143, 1)          1.0000
    (146, 1)          1.0000
    (147, 1)          1.0000
    (148, 1)          1.0000
    (149, 1)          1.0000
    (152, 1)          1.0000
    (153, 1)          1.0000
    (154, 1)          1.0000
    (155, 1)          1.0000
    (156, 1)          1.0000
    (158, 1)          1.0000
    (159, 1)          1.0000
    (160, 1)          1.0000
    (164, 1)          1.0000
    (176, 1)          1.0000
    (177, 1)          1.0000
    (178, 1)          1.0000
    (179, 1)          1.0000
    (180, 1)          1.0000
    (185, 1)          1.0000
    (186, 1)          1.0000
    (105, 2)          1.0000
    (106, 2)          1.0000
    (107, 2)          1.0000
    (108, 2)          1.0000
    (109, 2)          1.0000
    (110, 2)          1.0000
    (111, 2)          1.0000
    (112, 2)          1.0000
    (129, 2)          1.0000
    (130, 2)          1.0000
    (131, 2)          1.0000
    (165, 2)          1.0000
    (166, 2)          1.0000
    (167, 2)          1.0000
    (168, 2)          1.0000
    (169, 2)          1.0000
    (170, 2)          1.0000
    (173, 2)          1.0000
    (175, 2)          1.0000
    (182, 2)          1.0000
    (183, 2)          1.0000
    (201, 2)          1.0000
      (5, 3)          1.0000
     (10, 3)          1.0000
     (12, 3)          1.0000
     (13, 3)          1.0000
     (20, 3)          1.0000
     (30, 3)          1.0000
     (32, 3)          1.0000
     (44, 3)          1.0000
     (46, 3)          1.0000
     (51, 3)          1.0000
     (55, 3)          1.0000
     (56, 3)          1.0000
     (66, 3)          1.0000
     (67, 3)          1.0000
     (68, 3)          1.0000
     (71, 3)          1.0000
     (78, 3)          1.0000
     (90, 3)          1.0000
     (92, 3)          1.0000
     (96, 3)          1.0000
     (97, 3)          1.0000
    (100, 3)          1.0000
    (102, 3)          1.0000
    (113, 3)          1.0000
    (118, 3)          1.0000
    (121, 3)          1.0000
    (125, 3)          1.0000
    (127, 3)          1.0000
    (128, 3)          1.0000
    (133, 3)          1.0000
    (135, 3)          1.0000
    (136, 3)          1.0000
    (141, 3)          1.0000
    (142, 3)          1.0000
    (144, 3)          1.0000
    (145, 3)          1.0000
    (150, 3)          1.0000
    (161, 3)          1.0000
    (162, 3)          1.0000
    (174, 3)          1.0000
      (0, 4)          1.0000
     (17, 4)          1.0000
     (18, 4)          1.0000
     (21, 4)          1.0000
     (36, 4)          1.0000
     (38, 4)          1.0000
     (39, 4)          1.0000
     (43, 4)          1.0000
     (47, 4)          1.0000
     (48, 4)          1.0000
     (50, 4)          1.0000
     (61, 4)          1.0000
     (62, 4)          1.0000
     (63, 4)          1.0000
     (64, 4)          1.0000
     (65, 4)          1.0000
     (69, 4)          1.0000
    (103, 4)          1.0000
    (104, 4)          1.0000
    (151, 4)          1.0000
    (157, 4)          1.0000
    (187, 4)          1.0000
    (188, 4)          1.0000
    (189, 4)          1.0000
    (163, 5)          1.0000
    (171, 5)          1.0000
    (172, 5)          1.0000
```