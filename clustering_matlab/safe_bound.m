function LB = safe_bound(blk, At, C, b, y, Z2, Bt, y2, l)

  mu = 1.1;
  Aty = sdpnalAtyfun(blk, At, y);
  Znew = ops(C, '-', Aty);
  if ~isempty(Bt)
    Bty = sdpnalAtyfun(blk, Bt, y2);
    Znew = ops(Znew, '-', Bty);
  end
  if ~isempty(Z2)
     Znew = ops(Znew, '-', Z2); 
  end

  LB0 = b'*y + l'*y2; 
  pert = 0; 

  eigtmp = eig(full(Znew{1})); 
  idx = find(eigtmp < -1e-6); 
  % Xbar = mu * max(eig(full(X{1})));
  Xbar = mu; 
  numneg = length(idx); 
  
  if (numneg) 
     pert = pert + Xbar * sum(eigtmp(idx)); 
  end
      
  LB = LB0 + pert; 
  
end
