function [best_bound, best_Xopt, best_flag, best_ineq, cp_iter, cp_flag, ...
    avg_n_pair, avg_n_triangle, avg_n_clique, best_Bcell, best_l] = cutting_plane_loop(...
    blk, C, At, b, L, Z1, Z2, Bt, flag, bound, X, y, Bcell, l, ...
    k, original_n, original_trace, options, max_cp_iter, ...
    cp_tol, global_ub, bb_tol, eps_ineq, eps_active, max_ineq, ...
    max_pair_ineq, pair_perc, max_triangle_ineq, triangle_perc)

best_flag = flag;

%  0         - sdp is solved successfully
%  1         - sdp is not solved successfully
% -1, -2, -3 - sdp is partially solved successfully

best_bound = bound;
best_Xopt = X{1};
best_Bcell = Bcell;
best_l = l;
best_ineq = size(Bcell, 2);

avg_n_pair = 0;
avg_n_triangle = 0;
avg_n_clique = 0;

cp_iter = 0;
cp_flag = -1;  

% -3 - worst bound
% -2 - sdp is not solved successfully
% -1 - maximum number of iterations
%  0 - no violated inequalities
%  1 - maximum number of inequalities
%  2 - node must be pruned
%  3 - cutting-plane tolerance

gap = (global_ub - (bound + original_trace)) / global_ub;
    
if gap <= bb_tol
    cp_flag = 2;
    return;
end

% cutting plane iterations
for i=1:max_cp_iter
    
    % remove inactive constraints
    [Bcell, l] = remove_inactive_constraints(Bcell, l, blk, Bt, best_Xopt, eps_active);
    
    % separate constraints
    [B_pair, n_pair] = separate_pair(best_Xopt, eps_ineq, max_pair_ineq, pair_perc);
    [B_triangle, n_triangle] = separate_triangle(best_Xopt, eps_ineq, max_triangle_ineq, triangle_perc);
    [B_clique, n_clique] = separate_clique(best_Xopt, k, eps_ineq, original_n);
    n_ineq = n_pair + n_triangle + n_clique;
        
    if n_ineq <= 1
        cp_flag = 0;
        break;
    end

    temp_Bcell = [B_pair, B_triangle, B_clique];
    Bcell = [Bcell, temp_Bcell];
    l_pair_triangle = zeros(n_pair + n_triangle, 1);
    l_clique = (1/(original_n - k + 1)) * ones(n_clique, 1);
    temp_l = [l_pair_triangle; l_clique];
    l = [l; temp_l];
    current_ineq = size(Bcell, 2);
    
    clear B_pair B_triangle B_clique temp_Bcell temp_l

    Bt = svec(blk, Bcell, 1);
    
    [~, X, ~, y, Z1, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options, ... 
        X, [], y, Z1, Z2, [], []);
    
    flag = info.termcode;
    bound = safe_bound(blk, At, C, b, y, Z2, Bt, y2, l); 
    
    cp_iter = cp_iter + 1;
       
    avg_n_pair = avg_n_pair + n_pair;
    avg_n_triangle = avg_n_triangle + n_triangle;
    avg_n_clique = avg_n_clique + n_clique;
    
    % sdp is not solved successfully
    if flag == 1
        cp_flag = -2;
        break;
    end
    
    gap = (global_ub - (bound + original_trace)) / global_ub;
    
    % disp(' GAP: ')
    % disp(gap)
    
    if gap <= bb_tol
        cp_flag = 2;
        % update cutting plane    
        best_flag = flag;
        best_bound = bound;
        best_Xopt = X{1};
        best_Bcell = Bcell;
        best_l = l;
        best_ineq = current_ineq;
        break;
    end
    
    % bound does not improve in the current iteration
    if bound <= best_bound
        cp_flag = -3;
        break;
    end
    
    cp_stop = abs((best_bound - bound) / (bound + original_trace));
    
    % cutting plane stop tolerance
    if cp_stop <= cp_tol
        cp_flag = 3;
        % update cutting plane    
        best_flag = flag;
        best_bound = bound;
        best_Xopt = X{1};
        best_Bcell = Bcell;
        best_l = l;
        best_ineq = current_ineq;
        break;
    end
    
    % update cutting plane    
    best_flag = flag;
    best_bound = bound;
    best_Xopt = X{1};
    best_Bcell = Bcell;
    best_l = l;
    best_ineq = current_ineq;
    
    % max cuts allowed
    if best_ineq >= max_ineq
        cp_flag = 1;
        break;
    end
        
end

avg_n_pair = avg_n_pair / (cp_iter + 1e-8);
avg_n_triangle = avg_n_triangle / (cp_iter + 1e-8);
avg_n_clique = avg_n_clique / (cp_iter + 1e-8);

end
