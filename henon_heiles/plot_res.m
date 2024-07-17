function plot_res(R)
    Lw = 1;
    figure()
    semilogy(R(1,:),'linewidth',Lw)
    hold on;
    %title('Residual norm')
    xlabel('iteration','interpreter','latex')
    ylabel('Residual norm','Interpreter','latex')
    set(gca,'fontsize',20)
    set(gcf,'color','w');
    grid on
    %ylim([1e-15,1e2])
    
    for j = 2:size(R,1)
        %pause(1)
        semilogy(R(j,:),'linewidth',Lw)
    end
end