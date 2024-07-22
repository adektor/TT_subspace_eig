
% ---------------------------------------------------------------- %
% Approximate TT subspace for Henon-Heiles Hamiltonian
% ---------------------------------------------------------------- %

clear all;

n = 28; % num. Hermite points per dimension
d = 3;  % dimension

Htt = hh_tt(n,d);
Htt = round(Htt,1e-10);         % compress MPO ranks

% estimate largest eig. & shift
v0 = tt_rand(n,d,1); k = 5;
lam_max = upper_eig(v0,Htt,k,1e-6,512);
Htt = Htt - abs(lam_max)*tt_eye(n,d);

% Subspace iteration parameters
k = 5;           % subspace dimension
maxiter = 100;   % maximum # of iterations

% Chebyshev filter
m = 6;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

% inital basis
V0 = random_TT_basis(n,d,32,k);

%% Low-rank subspace iteration
tol = 1e-12;     % truncation tolerance
rmax = 4;        % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub,lamsub,R,cpu_t,~,~,Y,RQ,gradRQ,PgradRQ] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m,1);

%% Full subspace iteration
%H = full(Htt);
%V0f = zeros(2^L,k); for i = 1:k; V0f(:,i) = full(V0{i}); end % TT --> full
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1);
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1,a,b,m);

%% Plots
plot_res(R)
plot_rank(Vsub);
plot_eigs(lamsub);
plot_cpu_t(cpu_t);
plot_ritz_coeffs(Y,2);
plot_RQ(RQ,gradRQ,PgradRQ)
