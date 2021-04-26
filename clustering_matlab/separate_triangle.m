function [B_triangle, n_ineq] = separate_triangle(X, eps, max_triangle_ineq, triangle_perc)
    
    rng(1727);
    
    n = size(X, 1);
    n_idx = randperm(n);
    
    B_triangle = cell(1, max_triangle_ineq);
    counter_triangle = 1;
    
    violations = zeros(max_triangle_ineq, 1);
    
    stop = false;
    
    % TRIANGLE INEQUALITIES
    
    for i=n_idx
        for j=i+1:n
            for t=j+1:n
                 if (i ~= t)
                     viol1 = 0.5*X(i,j) + 0.5*X(j,i) + 0.5*X(i,t) + 0.5*X(t,i) - X(i,i) - 0.5*X(j,t) - 0.5*X(t,j);
                     viol2 = 0.5*X(i,j) + 0.5*X(j,i) + 0.5*X(j,t) + 0.5*X(t,j) - X(j,j) - 0.5*X(i,t) - 0.5*X(t,i);
                     viol3 = 0.5*X(i,t) + 0.5*X(t,i) + 0.5*X(j,t) + 0.5*X(t,j) - X(t,t) - 0.5*X(i,j) - 0.5*X(j,i);
                     if viol1 >= eps
                        violations(counter_triangle) = viol1;
                        B_triangle{counter_triangle} = sparse([i, j, i, t, i, j, t], [j, i, t, i, i, t, j], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n, n);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     if viol2 >= eps
                        violations(counter_triangle) = viol2;
                        B_triangle{counter_triangle} = sparse([i, j, j, t, j, i, t], [j, i, t, j, j, t, i], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n, n);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     if viol3 >= eps
                        violations(counter_triangle) = viol3;
                        B_triangle{counter_triangle} = sparse([i, t, j, t, t, i, j], [t, i, t, j, t, j, i], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n, n);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     
                     if counter_triangle - 1 == max_triangle_ineq 
                         stop = true;
                         break;
                     end   
                 end
            end
            if stop == true
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    n_ineq = counter_triangle - 1;
    max_ineq = max_triangle_ineq * triangle_perc;
    if n_ineq <= max_ineq
        B_triangle = B_triangle(1:n_ineq);
    else    
        violations = violations(1:n_ineq);
        [~, id_sorted] = sort(violations, 'descend');
        B_triangle = B_triangle(id_sorted);
        n_ineq = floor(max_ineq);
        B_triangle = B_triangle(1:n_ineq);
        
        clear id_sorted
    end
    
    clear n_idx violations 

end