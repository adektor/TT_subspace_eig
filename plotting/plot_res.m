function plot_res(R)

% Plot residual norm versus iteration

    Lw = 1;
    figure()
    semilogy(R(1,:),'linewidth',Lw)
    hold on;
    xlabel('iteration','interpreter','latex')
    ylabel('Residual norm','Interpreter','latex')
    set(gca,'fontsize',20)
    set(gcf,'color','w');
    grid on
    
    for j = 2:size(R,1)
        %pause(1)
        semilogy(R(j,:),'linewidth',Lw)
    end
end