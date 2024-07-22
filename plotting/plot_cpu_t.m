function plot_cpu_t(cpu_t)

% Plot the CPU-time of each iteration
Lw = 1.5;
figure()
semilogy(cpu_t,'linewidth',Lw)
xlabel('iteration','interpreter','latex')
ylabel('CPU time (s)','Interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on