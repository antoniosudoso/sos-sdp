function [Bcell, l] = remove_inactive_constraints(Bcell, l, blk, Bt, Xopt, tol)

    if isempty(Bcell) && isempty(l)
        return;
    end

    Xvec = svec(blk, Xopt, 1);
    active_id = abs(Bt{1}' * Xvec) <= tol;
    Bcell = Bcell(active_id);
    l = l(active_id);

end
