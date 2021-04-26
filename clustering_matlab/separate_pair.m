function [B_pair, n_ineq] = separate_pair(X, eps, max_pair_ineq, pair_perc)

    rng(1727);
    
    n = size(X, 1);
    n_idx = randperm(n);
    
    B_pair = cell(1, max_pair_ineq);
    counter_pair = 1;
    
    violations = zeros(max_pair_ineq, 1);
    
    stop = false;
    
    % PAIR INEQUALITIES
    
    for i=n_idx
        for j=i+1:n
            viol1 = 0.5*X(i, j) + 0.5*X(j, i) - X(i, i);
            viol2 = 0.5*X(i, j) + 0.5*X(j, i) - X(j, j);
            if viol1 >= eps
                violations(counter_pair) = viol1;
                B_pair{counter_pair} = sparse([i ; j ; i], [j; i; i], [-0.5; -0.5; 1], n, n);
                counter_pair = counter_pair + 1;
            end
            if viol2 >= eps
                violations(counter_pair) = viol2;
                B_pair{counter_pair} = sparse([i ; j ; j], [j; i; j], [-0.5; -0.5; 1], n, n);
                counter_pair = counter_pair + 1;
            end
            if counter_pair - 1 == max_pair_ineq
                stop = true;
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    n_ineq = counter_pair - 1;
    max_ineq = max_pair_ineq * pair_perc;
    if n_ineq <= max_ineq
        B_pair = B_pair(1:n_ineq);
    else
        violations = violations(1:n_ineq);
        [~, id_sorted] = sort(violations, 'descend');
        B_pair = B_pair(id_sorted);
        n_ineq = floor(max_ineq);
        B_pair = B_pair(1:n_ineq);
        
        clear id_sorted
    end
    
    clear violations
    
end
