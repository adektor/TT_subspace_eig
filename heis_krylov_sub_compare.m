

% -------------------------------- %
% Compare convergence of low-rank TT 
%   Lanczos and subspace methods 
% -------------------------------- %

clear all;

L = 10;   % # spins
J = 1;   % interaction strength

Htt = heis_tt(L,J);             % Hamiltonian MPO form
Htt = round(Htt,1e-12);         % compress MPO ranks

tol = 1e-12;     % truncation tolerance
rmax = 6;       % max rank

%% estimate largest eig. & shift
k = 8;
v0 = tt_rand(2,L,rmax); v0 = round(v0,tol,rmax); 
lam_max = upper_eig(v0,Htt,k,1e-6,rmax);
Htt = Htt - abs(lam_max)*tt_eye(2,L);

% Exact eigenvalues (ensure L small enough)
[lam,psi] = exact_eigs(Htt);

%% ------------------- Subspace ------------------- %% 
k = 10;           % subspace dimension
maxiter = 50;   % maximum # of iterations

% Chebyshev filter parameters
m = 2;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(2,L,rmax,k);
[Vsub,lam_sub,R,cpu_sub,a,b] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

% Subspace eigenvalue error
neigs = 5;
Es = zeros(neigs,maxiter);
for i = 1:neigs
    for j = 1:maxiter
        Es(i,j) = min(abs( lam(i) - lam_sub(:,j)));
    end
end
%%
figure()
semilogy(Es(1,:),'linewidth',Lw)
hold on;
%title(sprintf('Subspace eigenvalue error rmax = %i',rmax))
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on
%
for j = 2:neigs
    pause(1)
    semilogy(Es(j,:),'linewidth',Lw)
end

%% ------------------- Lanczos ------------------- %% 
maxiter = 100; % number of Krylov vectors
v0 = tt_rand(2,L,rmax); v0 = round(v0,tol,rmax); 

%[V_lan,T,eigs_lan] = lanczos_lr(Htt,v0,maxiter,tol,rmax,a,b,m);
[V_lan,T,eigs_lan] = lanczos_lr(Htt,v0,maxiter,tol,rmax);


% Krylov eigenvalue error
neigs = 5;
Ek = zeros(neigs,maxiter);
for i = 1:neigs
    for j = 2:maxiter-1
        Ek(i,j) = min(abs( lam(i) - eigs_lan{j}(:)));
    end
end

figure()
semilogy(Ek(1,:),'linewidth',Lw)
hold on;
%title(sprintf('Lanczos eigenvalue error rmax = %i',rmax))
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

for j = 2:neigs
    pause(1)
    semilogy(Ek(j,:),'linewidth',Lw)
end