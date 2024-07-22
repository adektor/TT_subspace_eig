function plot_RQ(RQ,gradRQ,PgradRQ)

% Plot Rayleigh quotient and gradient info vs. iteration

maxiter = length(RQ);
ngradRQ = zeros(1,maxiter);
nPgradRQ = zeros(1,maxiter);

for i = 1:maxiter
    ngradRQ(i) = norm(gradRQ{i});
    nPgradRQ(i) = norm(PgradRQ{i});
end

Lw = 1;

figure() % Gradients
semilogy(ngradRQ,'linewidth',Lw)
hold on;
semilogy(nPgradRQ,'--','linewidth',Lw)
xlabel('iteration','interpreter','latex')
ylabel('RQ gradient norm','Interpreter','latex')
legend('Euclidean', 'Riemannian','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

