function plot_eig_err(lam_approx,lam)

% Plot error between approximate eigenvalues and 
% benchmark eigenvalues

[k,maxiter] = size(lam_approx);

figure()
Lw = 1.5;
semilogy(min(abs(lam(:)-real(lam_approx(1,:)))),'linewidth',Lw)
hold on;
title('Eigenvalue error')
xlabel('iteration','interpreter','latex')
ylabel('error','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

for j = 2:k
    %pause(1)
    semilogy(min(abs(abs(lam(:))-abs(lam_approx(j,:)))),'linewidth',Lw)
end
