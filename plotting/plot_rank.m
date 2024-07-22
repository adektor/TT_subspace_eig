function plot_rank(Vsub)

% Plot the rank of subspace vectors versus iteration
% 
% Vsub{j}{i} TT vector i from iteration j
%       OR 
% Vsub{j}(:,i) full vector i from iteration j

k = length(Vsub{1});
maxiter = length(Vsub);
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
Lw = 1.5;
plot(r(1,:),'linewidth',Lw)
hold on
xlabel('iteration','interpreter','latex')
ylabel('$\left\| {\bf r} \right\|_1$','interpreter','latex')
set(gca,'fontsize',20)
set(gcf,'color','w');
grid on
for i = 2:k
    plot(r(i,:),'linewidth',Lw)
    %pause(1)
end