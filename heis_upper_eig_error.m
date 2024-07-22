
clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
Lw = 1.5; Ms = 10;

L = 10;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

tol = 1e-12;     % truncation tolerance
rmax = 32;       % max rank

%% estimate largest eig. & shift
k = 4; % # steps of Krylov
v0 = tt_rand(2,L,rmax); v0 = round(v0,tol,rmax); 
lam_max = upper_eig(v0,Htt,k,1e-6,rmax);

% Exact eigenvalues (ensure L small enough)
[lam,psi] = exact_eigs(Htt);

[l,i] = max(abs(lam));
fprintf('Largest eig (absolute value): %.3f, low-rank upper bound %.3f \n',l,lam_max)
