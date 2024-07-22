% ---------------------------------------------------------------- %
% Compare convergence of low-rank TT Lanczos and subspace methods 
% with corresponding methods using no low-rank truncation 
% ---------------------------------------------------------------- %

clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
Lw = 1.5; Ms = 10;

L = 10;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

%% Krylov

maxiter = 100;   % number of Krylov vectors
tol = 1e-12;     % truncation tolerance

E_kry = zeros(4,maxiter);
ranks_kry = zeros(4,maxiter);

v0 = tt_rand(2,L,1); v0 = round(v0,tol,1);
[V,T,eigs] = lanczos(full(Htt),full(v0),maxiter);
%%
pp = 1;
%for rmax = [32 30 16 4]
    rmax = 16; 
    [V_lr,T,eigs_lr] = lanczos_lr(Htt,v0,maxiter,tol,rmax);
    V_lr_full = lr_to_full_basis(V_lr);
%
    for jj = 1:maxiter
        ranks_kry(pp,jj) = max(V_lr{jj}.r);
        E_kry(pp,jj) = subspace(V(:,1:jj),V_lr_full(:,1:jj));
    end
%    pp = pp+1;
%end

% Subspace iteration

%maxiter = 100;   % # of subspace iterations
tol = 1e-10;     % truncation tolerance
k = 1;           % subspace dimension

E_sub = zeros(4,maxiter);
ranks_sub = zeros(4,maxiter);
pp = 1;

v0 = random_TT_basis(2,L,1,k);
v0f = zeros(2^L,k); for i = 1:k; v0f(:,i) = full(v0{i}); end % TT --> full
[V,~,~] = subspace_iter(v0f,full(Htt),maxiter,1);

%for rmax = [32 30 16 4]
    [V_lr,~,~,~] = subspace_iter_lr(v0,Htt,maxiter,tol,rmax);
    for jj = 1:maxiter
        ranks_sub(pp,jj) = max(V_lr{jj}{1}.r);
        E_sub(pp,jj) = subspace(V{jj},lr_to_full_basis(V_lr{jj}));
    end
%    pp = pp+1;
%end

%% Plots
LW = 1.3;
figure()
semilogy(E_sub(1,:),'linewidth',LW)
hold on
semilogy(E_kry(1,:),'--','linewidth',LW)

xlabel('iteration','interpreter','latex')
ylabel('$\angle(V,V_{\bf r})$','Interpreter','latex')
%title('angle between subspaces','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on
legend('Power iteration','Krylov','interpreter','latex','NumColumns',2)

ylim([1e-16,5e0])

figure()
plot(ranks_sub(1,:),'linewidth',LW)
hold on
plot(ranks_kry(1,:),'--','linewidth',LW)


xlabel('iteration','interpreter','latex')
ylabel('rank','interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on
%legend('$r_{\mathrm{max}} = 32$','$r_{\mathrm{max}} = 30$','$r_{\mathrm{max}} = 16$','$r_{\mathrm{max}} = 4$')

ylim([0,20])