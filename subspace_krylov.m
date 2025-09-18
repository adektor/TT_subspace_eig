% -------------------------------- %
% Compare convergence of low-rank TT 
%   Lanczos and subspace iteration  
% -------------------------------- %

clear all;
setup;

L = 10;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

tol = 1e-12;     % truncation tolerance
rmax = 6;        % max rank

%% estimate largest eig. & shift
k = 8;
v0 = tt_rand(2,L,rmax); v0 = round(v0,tol,rmax); 
lam_max = upper_eig(v0,Htt,k,1e-6,rmax);
Htt = Htt - abs(lam_max)*tt_eye(2,L);

% Exact eigenvalues (ensure L small enough)
[lam,psi] = exact_eigs(Htt);

%% ------------------- Subspace ------------------- %% 
k = 5;           % subspace dimension
maxiter = 400;   % maximum # of iterations

% Chebyshev filter parameters
m = 2;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(2,L,rmax,k);
[Vsub,lam_sub,R,cpu_sub,a,b] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m,0,0,1);

% Subspace eigenvalue error
neigs = 5;
Es = zeros(neigs,maxiter);
for i = 1:neigs
    for j = 1:maxiter
        Es(i,j) = min(abs( lam(i) - lam_sub(:,j)));
    end
end

%% Plot subspace error
Lw = 1.5;
figure()
semilogy(Es(1,:),'linewidth',Lw)
hold on;
%title(sprintf('Subspace eigenvalue error rmax = %i',rmax))
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',16)
set(gcf,'color','w');
grid on
ylim([1e-15 1e2])
%
for j = 2:neigs
    semilogy(Es(j,:),'linewidth',Lw)
end

%% ------------------- Lanczos ------------------- %% 
maxiter = 400; % number of Krylov vectors
v0 = tt_rand(2,L,rmax); v0 = round(v0,tol,rmax); 
[V_lan,T,eigs_lan] = lanczos_lr(Htt,v0,maxiter,tol,rmax);

% Krylov eigenvalue error
neigs = 5;
Ek = zeros(neigs,maxiter);
for i = 1:neigs
    for j = 2:maxiter-1
        Ek(i,j) = min(abs(lam(i) - eigs_lan{j}(:)));
    end
end

%% plot Krylov error
figure()
semilogy(Ek(1,:),'linewidth',Lw)
hold on;
%title(sprintf('Lanczos eigenvalue error rmax = %i',rmax))
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',16)
set(gcf,'color','w');
ylim([1e-15 1e2])
grid on

for j = 2:neigs
    semilogy(Ek(j,:),'linewidth',Lw)
end


%% Rank of exact eigenvectors and krylov vectors
[V,T,eigs] = lanczos(full(Htt),full(v0),100);

krylov_ranks = [];
for i = 1:20
    v_i = reshape(V(:,i),v0.n');
    v_i = tt_tensor(v_i);
    v_i = round(v_i,1e-10);
    krylov_ranks = [krylov_ranks, max(v_i.r)];
end

psi_ranks=[];
for i = 1:20
    psi_i = reshape(psi(:,i),v0.n');
    psi_i = tt_tensor(psi_i);
    psi_i = round(psi_i,1e-10);
    psi_ranks = [psi_ranks, max(psi_i.r)];
end

Lw = 1.5;
figure()
plot(psi_ranks,'o','linewidth',Lw)
hold on;
plot(krylov_ranks,'x','linewidth',Lw)
xlabel('j','interpreter','latex')
ylabel('$\mathrm{max}({\bf r})$','Interpreter','latex')
set(gca,'fontsize',16)
set(gcf,'color','w');
grid on
legend('$\psi_j$ eigenvectors','$v_j$ Krylov vectors','interpreter','latex')
ylim([0,40])
