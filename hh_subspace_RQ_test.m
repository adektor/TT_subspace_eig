
% ---------------------------------------------------------------- %
% Approximate TT subspace for Henon-Heiles Hamiltonian
% A script for checking Rayleigh quotient for various ranks. 
% ---------------------------------------------------------------- %

clear all;

n = 56; % num. Hermite points per dimension
d = 2;  % dimension

Htt = hh_tt(n,d);
Htt = round(Htt,1e-10);         % compress MPO ranks

% estimate largest eig. & shift
v0 = tt_rand(n,d,1); k = 5;
lam_max = upper_eig(v0,Htt,k,1e-6,512);
Htt = Htt - abs(lam_max)*tt_eye(n,d);

[lam,psi] = exact_eigs(Htt); lam = real(lam);
psi1 = reshape(psi(:,1),n,n); psi1tt = tt_tensor(psi1);

% Subspace iteration parameters
k = 3;           % subspace dimension
maxiter = 100;   % maximum # of iterations

% Chebyshev filter
m = 7;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

% inital basis
V0 = random_TT_basis(n,d,32,k);

%% Low-rank subspace iteration
tol = 1e-12;     % truncation tolerance
RQrsub = [];
RQr = [];
bnd = [];
ev_er1 = [];
ev_er2 = [];

ranks = 1:1:10;
for rmax = ranks
    for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
    [Vsub,lamsub,R,cpu_t,~,~,Y,RQ,gradRQ,PgradRQ] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m,1);
    RQrsub = [RQrsub,real(RQ(end))];

    psi1_r = round(psi1tt,tol,rmax);
    psi1_r = psi1_r./norm(psi1_r); 
    RQr = [RQr, dot(psi1_r,Htt*psi1_r)];
    
    eps_r = norm(psi1_r-psi1tt);
    bnd = [bnd, lam(end) + (lam(1)-lam(end))*(1-eps_r^2)];
    
    ev_er1 = [ev_er1, subspace(full(psi1tt),full(Vsub{end}{1}))];
    ev_er2 = [ev_er2, subspace(full(psi1tt),full(psi1_r))];
end
%%
figure() % eigenvalue error
semilogy(ranks,abs(RQrsub-lam(1)),'LineWidth',1.5);
hold on
semilogy(ranks,abs(RQr-lam(1)),'--','LineWidth',1.5);
semilogy(ranks,abs(bnd-lam(1)),'-','color', [.5 .5 .5]);
%semilogy(ranks,lam(1)*ones(length(ranks),1),'color','k','linewidth',1.5);

xlabel('rank','interpreter','latex');
ylabel('$|\rho(x)-\lambda_1|$','interpreter','latex')
legend('$\rho(\psi_r)$','$\rho(T_r(\psi))$','upper bound')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

figure() % eigenvector error
semilogy(ranks,ev_er1,'LineWidth',1.5);
hold on
semilogy(ranks,ev_er2,'--','LineWidth',1.5);

xlabel('rank','interpreter','latex');
%ylabel('$\angle (\psi-\psi_r)$','interpreter','latex')
legend('$\angle (\psi-\psi_r)$','$\angle (\psi-T_r(\psi))$')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

%% Full subspace iteration
%H = full(Htt);
%V0f = zeros(2^L,k); for i = 1:k; V0f(:,i) = full(V0{i}); end % TT --> full
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1);
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1,a,b,m);

%% Plots
 plot_res(R)
% plot_rank(Vsub);
% plot_eigs(lamsub);
% plot_cpu_t(cpu_t);
% plot_ritz_coeffs(Y,2);
% animate_ritz_coeffs(Y,2);
% plot_RQ(RQ,gradRQ,PgradRQ)

