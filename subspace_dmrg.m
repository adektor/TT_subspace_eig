% -------------------------------- %
% Compare convergence of low-rank TT 
%   subspace iteration and DMRG 
% -------------------------------- %

clear all;

L = 100;                        % # spins
n = 3;                          % mode sizes
Htt = heis_tt_spin1(L);         % Hamiltonian MPO form
Htt = round(Htt,1e-10);         % compress MPO ranks

% estimate largest eig. & shift
v0 = tt_rand(n,L,1); m = 5;
lam_max = upper_eig(v0,Htt,m,1e-6,512);
Htt = Htt - abs(lam_max)*tt_eye(n,L);

%% DMRG
rmax = 100;         % max rank
tolDMRG = 1e-12;    % truncation tolerance
n_eigvecs = 1;      % # of eigenvectors

[x_dmrg, E_dmrg] = dmrg_eig(Htt, tolDMRG, ...
                        'b', n_eigvecs, ...
                        'nswp', 8, ...
                        'max_full_size', 1024, ...
                        'rmax', rmax, ...
                        'verb',1);

%% Subspace
% Subspace iteration parameters
m = 2;           % subspace dimension
maxiter = 500;   % maximum # of iterations

% Chebyshev filter parameters
k = 4;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

% truncation parameters
rmax = 100;      % max rank
tol = 1e-12;     % tolerance

V0 = random_TT_basis(n,L,rmax,m);
[Vsub,lamsub,R,cpu_t] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,k,1,lam_max,1);


%% Print energies
E_ref = -140.14840390392;
fprintf("reference energy: %.6f \n", E_ref)
fprintf("DMRG energy: %.6f \n", E_dmrg + lam_max)
fprintf("subspace energy: %.6f \n", lamsub(end) + lam_max)

