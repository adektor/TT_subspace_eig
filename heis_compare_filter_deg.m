clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
Lw = 1.5; Ms = 10;

L = 16;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

%Httf = full(Htt);              % full Hamiltonian for validation.
%H = heis(L,J);     
%norm(H-Httf,'fro')

% estimate largest eig. & shift
v0 = tt_rand(2,L,1); k = 5;
lam_max = upper_eig(v0,Htt,k,1e-6,16);
Htt = Htt - abs(lam_max)*tt_eye(2,L);

% Subspace iteration parameters
k = 10;           % subspace dimension
maxiter = 50;   % maximum # of iterations

V0 = random_TT_basis(2,L,512,k);

tol = 1e-12;     % truncation tolerance
rmax = 8;        % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end

%% No filter
[Vsub0,~,R0,cpu_t0] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax);

%% Degree 2 
m = 2;                 % degree
a = -1; b = 0;         % window of spectrum to avoid
[Vsub2,~,R2,cpu_t2] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Degree 8 
m = 8;                 % degree
a = -1; b = 0;         % window of spectrum to avoid
[Vsub8,~,R8,cpu_t8] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Degree 16
m = 16;                 % degree
a = -1; b = 0;         % window of spectrum to avoid
[Vsub16,~,R16,cpu_t16] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);
