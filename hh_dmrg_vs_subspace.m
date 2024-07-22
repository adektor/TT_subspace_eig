
% ---------------------------------------------------------------- %
% Compare block-DMRG and approximate TT subspace for Henon-Heiles
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

%% Block DMRG
% parameters
tolDMRG = 1e-5;              % truncation tolerance
n_eigvecs = 3;           % # of eigenvectors

[x,theta,testdata] = dmrg_eig(Htt, tolDMRG, 'b', n_eigvecs, 'max_full_size', 1024);

% move block around and compute residuals/memory:
[R_dmrg,numel_blk] = check_block_dmrg(x,theta,Htt,tolDMRG*1e-1);

%% Subspace iteration
k = 8;           % subspace dimension
maxiter = 40;   % maximum # of iterations
rmax = 20;       % max rank

% Chebyshev filter parameters
m = 6;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(n,d,8,k);

tol = 1e-12;     % truncation tolerance
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub,lamsub,R,cpu_t] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

R_sub = [];
for i = 1:k
    R_sub = [R_sub, norm(Htt*Vsub{end}{i} - lamsub(i,end)*Vsub{end}{i})];
end

% Count number of elements
num_el_sub=0;
for i = 1:n_eigvecs
    Vsub{end}{i} = round(Vsub{end}{i},tolDMRG);
    num_el_sub = num_el_sub + numel(Vsub{end}{i}.core);
end 

fprintf('Number of degrees of freedom to represent the first %i eigenvectors up to tolerance %.1e \n subspace: %.3f \t block-DMRG %.3f \n',n_eigvecs,tolDMRG,num_el_sub,min(numel_blk))