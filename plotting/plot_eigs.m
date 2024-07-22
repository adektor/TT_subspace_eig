function plot_eigs(lamsub,lam)

% Plot eigenvalues versus iteration
% 
% INPUTS: 
% lamsub(j,i) approximate eigenvalue j from iteration i
% lam(j) actual eigenvalue j

[k,maxiter] = size(lamsub);

figure()
Lw = 1.5;
if nargin == 2 % plot true eigenvalues
    for j = 1:k
        plot(ones(1,maxiter)*lam(j),"Color", [0.4, 0.4, 0.4, .5],'linewidth',Lw)
        hold on
        title('Eigenvalues')
    end
end

for j = 1:k
    plot(real(lamsub(j,:)),'linewidth',Lw)
    hold on
    %pause(1)
end
xlabel('iteration','interpreter','latex')
ylabel('$\lambda_i$','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on