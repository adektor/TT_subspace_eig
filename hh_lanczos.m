
% ---------------------------------------------------------------- %
% Approximate TT Lanczos for Henon-Heiles Hamiltonian
% ---------------------------------------------------------------- %

clear all;

n = 28; % num. Hermite points per dimension
d = 2;  % dimension

Htt = hh_tt(n,d);
Htt = round(Htt,1e-10);         % compress MPO ranks

rmax = 16;
tol = 1e-12;
v0 = tt_rand(n,d,1); v0 = round(v0,tol,rmax);
maxiter = 100;

% low-rank vectors
[V_lan,T,eigs_lan] = lanczos_lr(Htt,v0,maxiter,tol,rmax);

% full vectors
%[V_lanf,Tf,eigs_lanf] = lanczos(H,V0f(:,1),maxiter);

% exact eigs (only for small d)
[lam,psi] = exact_eigs(Htt); 

% extract first k eigenvalues at each iteration > k
k=10;
lam_lanczos = [];
for j = k:maxiter
    lam_lanczos = [lam_lanczos, eigs_lan{j}(1:k)];
end

%% Plot 
plot_eig_err(lam_lanczos,lam)
