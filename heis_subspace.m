
% -------------------------------- %
% Low-rank TT subspace iteration 
%  w/ polynomial acceleration 
%  for Heisenberg Hamiltonian 
% -------------------------------- %

clear all;

L = 12;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

% estimate largest eig. & shift
v0 = tt_rand(2,L,1); k = 5;
lam_max = upper_eig(v0,Htt,k,1e-6,512);
Htt = Htt - abs(lam_max)*tt_eye(2,L);

% Subspace iteration parameters
k = 5;           % subspace dimension
maxiter = 100;   % maximum # of iterations

% Chebyshev filter parameters
m = 5;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(2,L,512,k);

tol = 1e-12;     % truncation tolerance
rmax = 2;        % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub,lamsub,R,cpu_t,~,~,Y,RQ,gradRQ,PgradRQ] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m,1);

%% full subspace iteration
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1);
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1,a,b,m);

%% Plot
plot_res(R)
plot_rank(Vsub);
plot_eigs(lamsub);
plot_cpu_t(cpu_t);
plot_ritz_coeffs(Y,2);
plot_RQ(RQ,gradRQ,PgradRQ)
