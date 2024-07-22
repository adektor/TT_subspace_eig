function animate_ritz_coeffs(Y,j)

% Plot the coefficients for the jth Ritz vector 
% versus iteration

[maxiter,k,~] = size(Y);

figure()
vid = VideoWriter('coefficients.avi');
open(vid);

for i = 1:maxiter
    semilogy(abs(Y(i,:,j)),'o','MarkerSize',8,'linewidth',1.5);
    ylim([1e-17,2e0])
    title(sprintf('iteration %i',i))
    xlabel('$\varphi_j$', 'Interpreter','Latex')
    ylabel('$j$', 'Interpreter','Latex')
    set(gca,'fontsize',20)
    set(gcf,'color','w');
    grid on
    frame = getframe(gcf);
    writeVideo(vid,frame);
    pause(0.075)
end
close(vid);

end