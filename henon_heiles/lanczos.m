function [V,T,eigs] = lanczos(A,v0,maxiter)

% Lanczos method 
% INPUT: A --> n x n matrix
%        v0 --> n x 1 vector with norm 1
%        maxiter --> number of Lanczos steps
%
% OUTPUT: V --> n x maxiter matrix containing Lanczos vectors
%         T --> tridiagonal matrix of Lanczos coefficients alpha,beta
%         eigs --> cell array containing approximate eigenvalues at each iteration

% initialize:
n = length(v0);
alpha = [];
beta = 0;

V = zeros(n,maxiter); V(:,1) = v0;
eigs = cell(1,maxiter);

% first vector obtained with two-term recurrence:
V(:,2) = A*V(:,1);
alpha(1) = V(:,2)'*V(:,1);
V(:,2) = V(:,2) - alpha(1)*V(:,1);

for j = 2:maxiter
    beta = [beta norm(V(:,j))];
    V(:,j) = V(:,j)./beta(j);
    
    % obtain the rest with three-term recurrence:
    V(:,j+1) = A*V(:,j);
    
    alpha = [alpha V(:,j+1)'*V(:,j)];
    V(:,j+1) = V(:,j+1) - alpha(j)*V(:,j) - beta(j)*V(:,j-1);

    % construct T and compute eigs
    T = diag(alpha) + diag(beta(2:end),-1) + diag(beta(2:end),1);
    eigs{j} = eig(T); eigs{j} = sort(eigs{j},'descend');

    if mod(j,10) == 0 
        fprintf("Iteration %i: \n",j)
    end
end

end