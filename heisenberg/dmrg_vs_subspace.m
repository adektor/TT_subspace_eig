clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
Lw = 1.5; Ms = 10;

num_el_dmrg = zeros(1,5);
num_el_sub = zeros(1,5);
pp = 1;
Llist = [64];
%for L = Llist
    %fprintf("%i spins \n",L)
L = 64;   % # spins
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

%% Block DMRG
% parameters
tol = 1e-9;
k = 15;           % # of eigenvectors 

[x,theta,testdata] = dmrg_eig(Htt, tol, 'b', k, 'max_full_size', 1024);

% move block around and compute residuals/memory:
[R_dmrg,numel_blk] = check_block_dmrg(x,theta,Htt,tol);

%% Subspace iteration
k = 15;           % subspace dimension
maxiter = 300;   % maximum # of iterations
rmax = 4;       % max rank

% Chebyshev filter parameters
m = 8;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(2,L,8,k);

tol = 1e-12;     % truncation tolerance
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub,lamsub,R,cpu_t] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%
R_sub = [];
for i = 1:k
    R_sub = [R_sub, norm(Htt*Vsub{end}{i} - lamsub(i,end)*Vsub{end}{i})];
end

% Count number of elements
for i = 1:k-1
    num_el_sub(pp) = num_el_sub(pp) + numel(Vsub{end}{i}.core);
end
pp = pp + 1;
%end

%% Plot
% figure()
% loglog(Llist,num_el_sub,'linewidth',1.4)
% hold on
% loglog(Llist,num_el_dmrg,'--','linewidth',1.4)
% xlim([Llist(1),Llist(end)])
% xticks(Llist)
% legend('subspace','block DMRG','location','northwest')
% xlabel('$L$','interpreter','latex')
% ylabel('num entries','interpreter','latex')
% %
% set(gca,'fontsize',20)
% set(gcf,'color','w');
% grid on