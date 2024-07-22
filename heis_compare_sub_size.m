clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
Lw = 1.5; Ms = 10;

L = 32;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

%Httf = full(Htt);              % full Hamiltonian for validation.
%H = heis(L,J);     
%norm(H-Httf,'fro')

% estimate largest eig. & shift
v0 = tt_rand(2,L,1); k = 5;
lam_max = upper_eig(v0,Htt,k,1e-6,512);
Htt = Htt - abs(lam_max)*tt_eye(2,L);

% Subspace iteration parameters
maxiter = 500;   % maximum # of iterations

m = 3;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(2,L,512,k);

tol = 1e-12;     % truncation tolerance
rmax = 2;        % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end

%% Dimension 2
k = 2;       
V0 = random_TT_basis(2,L,512,k);
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end

[Vsub2,~,R2,cpu_t2] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Dimension 4 
k = 4;       
V0 = random_TT_basis(2,L,512,k);
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end

[Vsub4,~,R4,cpu_t4] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Dimension 8 
k = 8;       
V0 = random_TT_basis(2,L,512,k);
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end

[Vsub8,~,R8,cpu_t8] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Dimension 16
k = 16;       
V0 = random_TT_basis(2,L,512,k);
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end

[Vsub16,~,R16,cpu_t16] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Plot
plot_res(R4)
plot_res(R8)
plot_res(R16)