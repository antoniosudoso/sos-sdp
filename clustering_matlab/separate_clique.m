function [B_clique, n_ineq] = separate_clique(X, k, eps, old_n)
    
    n = size(X, 1);
    
    B_clique = cell(1, n);
    counter_clique = 1;
    
    % CLIQUE INQUALITIES
    
    V = (1:n)';
    
    for j=1:n
        C = j;
        C_size = size(C, 1);
        while C_size <= k
            V_minus_C = setdiff(V, C);
            temp_sum = [];
            temp_v = [];
            for v=1:n
                sum = 0;
                if ismember(v, V_minus_C)
                    for i=1:C_size
                        sum = sum + X(C(i), v);
                    end
                    temp_sum = [temp_sum; sum];
                    temp_v = [temp_v; v];
                end    
            end
            [~, argmin] = min(temp_sum);
            u = temp_v(argmin);
            C = [C; u];  
            C_size = size(C, 1);    
        end 
        
        
        lhs = 0;
        row_idx = [];
        col_idx = [];
        for i=1:C_size
            for z=1:C_size
                if i < z
                    lhs = lhs + X(C(i), C(z));
                    row_idx = [row_idx; C(i)];
                    col_idx = [col_idx; C(z)];
                end
            end
        end
        
        size_idx = size(row_idx, 1);
        viol = lhs - (1 / (old_n - k + 1));
               
        if viol < -eps
            B_clique{counter_clique} = sparse([row_idx; col_idx], [col_idx; row_idx], 0.5 * ones(2 * size_idx, 1), n, n);
            counter_clique = counter_clique + 1;
        end
                    
    end
    
    n_ineq = counter_clique - 1;
    B_clique = B_clique(1:n_ineq);
    
    clear V
    
end
