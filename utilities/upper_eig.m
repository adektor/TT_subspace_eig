function lam_max = upper_eig(v0,H,k,tol,rmax)

% Estimate the largest eigenvalue using a few steps of Lanczos method. 
% INPUT: v0 --> either full vector or TT tensor
%        H --> Hamiltonian 
%        k --> # of iterations (small)
%        tol --> truncation tolerance
%        rmax --> maximum rank

% OUTPUT: lam_max --> approximate upper bound for dominant eigenvalue

if isa(v0,'double')
    [V_lan,T,eigs_lan] = lanczos(H,v0,k);
    lam_max = max(abs(eigs_lan{end})) + norm(V_lan(:,end));
else
    [V_lan,T,eigs_lan] = lanczos_lr(H,v0,k,tol,rmax);
    lam_max = max(abs(eigs_lan{end})) + norm(V_lan{end});
end

end