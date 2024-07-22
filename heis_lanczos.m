
% -------------------------------- %
% Low-rank TT approximate Lanczos 
%   for Heisenberg Hamiltonian 
% -------------------------------- %

clear all;

L = 6;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

rmax = 2;
tol = 1e-12;
v0 = tt_rand(2,L,1); v0 = round(v0,tol,rmax);
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
