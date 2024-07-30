
% ------------------------------------------------------- %
% Low-rank TT subspace iteration and DMRG
%  w/ polynomial acceleration 
%  for Heisenberg spin-1 Hamiltonian 
% ------------------------------------------------------- %

clear all;

L = 10;   % # spins
n = 3; 
Htt = heis_tt_spin1(L);             % Hamiltonian MPO form
Htt = round(Htt,1e-10);             % compress MPO ranks

%[lam,psi] = exact_eigs(Htt);

%%
% estimate largest eig. & shift
v0 = tt_rand(n,L,1); k = 5;
lam_max = upper_eig(v0,Htt,k,1e-6,512);
Htt = Htt - abs(lam_max)*tt_eye(n,L);

%% DMRG

% parameters
rmax = 100;        % max rank
tolDMRG = 1e-5;              % truncation tolerance
n_eigvecs = 1;               % # of eigenvectors

[x,theta,testdata] = dmrg_eig(Htt, tolDMRG, 'b', n_eigvecs,'nswp',50,'max_full_size', 1024,'rmax',rmax,'verb',1);
[R_dmrg,numel_blk] = check_block_dmrg(x,theta,Htt,tolDMRG*1e-3);

%% Subspace
% Subspace iteration parameters
k = 2;           % subspace dimension
maxiter = 30;   % maximum # of iterations

% Chebyshev filter parameters
m = 5;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

rmax = 100;      % max rank
tol = 1e-4;     % truncation tolerance

V0 = random_TT_basis(n,L,rmax,k);
[Vsub,lamsub,R,cpu_t,~,~,Y,RQ,gradRQ,PgradRQ] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,...
                                                        'a',a, ...
                                                        'b',b, ... 
                                                        'm',m, ...
                                                        'var_mv',1 );
plot_res(R)

%% Print residuals
fprintf("Residuals: \n block-DMRG: %.2e \n subspace: %.2e \n", min(R_dmrg),R(end,end))

%% full subspace iteration
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1);
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1,a,b,m);

%% Plot
% plot_res(R)
% plot_rank(Vsub);
% plot_eigs(lamsub);
% plot_cpu_t(cpu_t);
% plot_ritz_coeffs(Y,2);
% animate_ritz_coeffs(Y,1);
% plot_RQ(RQ,gradRQ,PgradRQ)
