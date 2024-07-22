function plot_ritz_coeffs(Y,j)

% Plot the coefficients for the jth Ritz vector 
% versus iteration

[maxiter,k,~] = size(Y);

figure()
for i = 1:k
    semilogy(abs(Y(:,i,j)));
    hold on;
end

hold on;
xlabel('iteration','interpreter','latex')
ylabel(sprintf('coefficients of Ritz vector %i',j),'Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on

end