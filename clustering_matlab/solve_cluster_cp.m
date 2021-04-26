function [best_bound, best_Xopt, best_flag, best_ineq, cp_iter, cp_flag, ...
    avg_n_pair, avg_n_triangle, avg_n_clique, best_Bcell, best_l] = solve_cluster_cp(n_threads, ...
    WWt, Acell, b, k, original_n, original_trace, ...
    tol, verbose, max_cp_iter, cp_tol, global_ub, bb_tol, ...
    eps_ineq, eps_active, max_ineq, max_pair_ineq, pair_perc, max_triangle_ineq, triangle_perc)

% note that Z1 is the psd matrix and Z2 is the non-negative matrix.

warning off;

% disp(' N_THREADS: ')
maxNumCompThreads(n_threads);
% disp(maxNumCompThreads)

n = size(WWt,1);

blk = cell(1);
blk{1, 1} = 's';
blk{1, 2} = n;

C = cell(1);
C{1} = -WWt;

At = svec(blk, Acell, 1);

L = cell(1);
L{1} = 0; % we want X >= 0

options.tol = tol;
options.printlevel = verbose;
% options.maxtime = 7200;
% options.stopoption = 0;

% solve without inequalities
[~, X, ~, y, Z1, Z2, ~, ~, info, ~] = admmplus(blk, At, C, b, L, [], [], [], [], options);
flag = info.termcode;
bound = safe_bound(blk, At, C, b, y, Z2, [], 0, 0);

Bcell = cell(1, 0);
l = zeros(0, 1);
Bt = svec(blk, Bcell, 1);

[best_bound, best_Xopt, best_flag, best_ineq, cp_iter, cp_flag, ...
    avg_n_pair, avg_n_triangle, avg_n_clique, best_Bcell, best_l] = cutting_plane_loop(...
    blk, C, At, b, L, Z1, Z2, Bt, flag, bound, X, y, Bcell, l, ...
    k, original_n, original_trace, options, max_cp_iter, cp_tol, ...
    global_ub, bb_tol, eps_ineq, eps_active, max_ineq, ...
    max_pair_ineq, pair_perc, max_triangle_ineq, triangle_perc);


clear C At Bcell Bt l X y Z1 Z2


end
