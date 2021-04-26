function [new_Bcell, new_l] = inherit_cuts(parent_Bcell, parent_l, n, parent_n, inherit_perc, i, j)
    
    if inherit_perc == 0
        new_Bcell = cell(1, 0);
        new_l = zeros(0, 1);
        return;
    end

    n_ineq = size(parent_Bcell, 2);
    
    max_inherit_ineq = floor(inherit_perc * n_ineq);
    
    if n == parent_n
        new_Bcell = parent_Bcell(1:max_inherit_ineq);
        new_l = parent_l(1:max_inherit_ineq);
        return;
    end
    
    new_Bcell = cell(1, max_inherit_ineq);
    new_l = zeros(max_inherit_ineq, 1);
    
    counter = 1;
    
    for c=1:n_ineq
        
        [id_i, id_j, v] = find(parent_Bcell{c});
        if any(id_i == i) || any(id_i == j) || any(id_j == i) || any(id_j == j)
            continue;
        end
        
        d = size(id_i, 1);
        
        for k=1:d
            
            if id_i(k) > j
                id_i(k) = id_i(k) - 1;
            end
            
            if id_j(k) > j
                id_j(k) = id_j(k) - 1;
            end
            
        end
        
        new_Bcell{counter} = sparse(id_i, id_j, v, n, n);
        new_l(counter) = parent_l(c);
        counter = counter + 1;
        
        if counter - 1 == max_inherit_ineq
            break;
        end
                
    end
    
    n_ineq = counter - 1;
    new_Bcell = new_Bcell(1:n_ineq);
    new_l = new_l(1:n_ineq);
    
end