clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
Lw = 1.5; Ms = 10;

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

% Subspace iteration parameters
k = 15;           % subspace dimension
maxiter = 1000;   % maximum # of iterations

% Chebyshev filter parameters
m = 8;                 % degree
a = -1; b = 0;         % window of spectrum to avoid

V0 = random_TT_basis(2,L,512,k);
%V0f = zeros(2^L,k); for i = 1:k; V0f(:,i) = full(V0{i}); end % TT --> full

% Rank 2
tol = 1e-12;     % truncation tolerance
rmax = 2;       % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub2,lamsub2,R2,cpu_t2] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

% Rank 32
tol = 1e-12;      % truncation tolerance
rmax = 32;       % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub32,lamsub32,R32,cpu_t32] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

% Rank 128
tol = 1e-12;      % truncation tolerance
rmax = 128;       % max rank
for i = 1:k; V0{i} = round(V0{i},tol,rmax); end
[Vsub128,lamsub128,R128,cpu_t128] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m);

%% Full iteration
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1);
%[Vsub,lamsub,R] = subspace_iter(V0f,H,maxiter,1,a,b,m);

%% Plot: ranks
r = zeros(k,maxiter);
for i = 1:k
    for j = 1:maxiter
        if isa(Vsub{j},'double')
            ten = reshape(Vsub{j}(:,i),ones(1,L)*2);
            tt = tt_tensor(ten,1e-11);
            r(i,j) = sum(tt.r);
        else
            r(i,j) = sum(Vsub{j}{i}.r);
        end
    end
end

figure()
plot(r(1,:),'linewidth',Lw)
hold on
%title('$\left\| {\bf r} \right\|_1$','interpreter','latex')
xlabel('iteration','interpreter','latex')
ylabel('$\left\| {\bf r} \right\|_1$','interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on
for i = 2:k
    plot(r(i,:),'linewidth',Lw)
    %pause(1)
end

%% Plot: eigenvalues versus iteration
figure()

% for j = 1:k
%     plot(ones(1,maxiter)*lam(j),"Color", [0.4, 0.4, 0.4, .5],'linewidth',Lw)
%     hold on
%     title('Eigenvalues')
% end

%ylim([-3.65 -3])
for j = 1:k
    plot(real(lamsub(j,:)),'linewidth',Lw)
    hold on
    %pause(1)
end
%ylim([-18.1 -17.75])
xlabel('iteration','interpreter','latex')
ylabel('$\lambda_i$','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

%% Plot: CPU-time of iterations
figure()
semilogy(cpu_t,'linewidth',Lw)
xlabel('iteration','interpreter','latex')
ylabel('CPU time (s)','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

%% Plot: eigenvalue error
figure()
semilogy(abs(lam(1)-real(lamsub(1,:))),'linewidth',Lw)
hold on;
title('Eigenvalue error')
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

for j = 2:k
    pause(1)
    semilogy(abs(lam(j)-real(lamsub(j,:))),'linewidth',Lw)
end

%% Plot: eigenspace error
E_eig = zeros(k,maxiter);
for j = 1:maxiter
    for i = 1:k-2
        if isa(Vsub{j},'double')
            E_eig(i,j) = subspace(Vsub{j}(:,i),psi(:,i));
        else
            E_eig(i,j) = subspace(full(Vsub{j}{i}),psi(:,i));
        end
    end
    if isa(Vsub{j},'double')
        E_eig(k-1,j) = subspace([Vsub{j}(:,9),Vsub{j}(:,10)],psi(:,9:10));
    else
        E_eig(k-1,j) = subspace([full(Vsub{j}{9}),full(Vsub{j}{10})],psi(:,9:10));
    end
end

figure()
semilogy(E_eig(1,:),'linewidth',Lw)
hold on;
title('Eigenspace error')
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on
ylim([1e-12,1e1])

for i = 2:size(E_eig,1)
    pause(1)
    semilogy(E_eig(i,:),'linewidth',Lw)
end



%% ------------------- Lanczos ------------------- %%
rmax = 32;
tol = 1e-12;
v0 = tt_rand(2,L,1); v0 = round(v0,tol,rmax);
maxiter = 100;

%%
%[V_lanf,Tf,eigs_lanf] = lanczos(H,V0f(:,1),maxiter);
[V_lan,T,eigs_lan,eigs_lan2] = lanczos_lr(Htt,V0{1},maxiter);

%% Extract first k eigenvalues versus iteration
k=10;
lam_lanczos = [];
for j = k:maxiter
    %lam_lanczos = [lam_lanczos, eigs_lan{j}(1:k)];
    lam_lanczos = [lam_lanczos, eigs_lanf{j}(1:k)];
end

%% Plot: eigenvalue error
figure()
semilogy(min(abs(lam(:)-real(lam_lanczos(1,:)))),'linewidth',Lw)
hold on;
title('Eigenvalue error')
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

for j = 2:k
    pause(1)
    semilogy(min(abs(abs(lam(:))-abs(lam_lanczos(j,:)))),'linewidth',Lw)
end

%% Plot: eigenvalues versus iteration
figure()
for j = 1:k
    plot(ones(1,maxiter-k+2)*real(lam(j)),"Color", [0.4, 0.4, 0.4, .5],'linewidth',Lw)
    hold on
    title('Eigenvalues')
end

xlabel('iteration','interpreter','latex')
ylabel('$\lambda_i$','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');

ylim([-4.2 -3.9])
for j = 1:k
    plot(-abs(lam_lanczos(j,:)),'linewidth',Lw)
    pause(.7)
end


%% Plot
[eigs2,acc2] = process_eigs(eigs_lan,lam);

% Lanczos eigs
figure()
pointsize = 50;
hold on
for i = 2:maxiter
    scatter(i*ones(1,i),eigs2{i},pointsize,acc2{i},'filled');
end
%set(gca,'yscale','log')
title('Lanczos')
xlabel('iteration')
ylabel('eigenvalues')
colorbar
%ylim([-12,12])
set(gca,'ColorScale','log')

set(gcf,'color','w')
set(gca,'fontsize',20)
grid on
box on

% Lanczos low-rank
figure()
pointsize = 50;
hold on
for i = 2:maxiter
    scatter(i*ones(1,i),real(eigs2_lr{i}),pointsize,acc2_lr{i},'filled');
end

%set(gca,'yscale','log')
title('Lanczos low-rank')
xlabel('iteration')
ylabel('eigenvalues')
colorbar
set(gca,'ColorScale','log')
%ylim([-12,12])

set(gcf,'color','w')
set(gca,'fontsize',20)
grid on
box on
