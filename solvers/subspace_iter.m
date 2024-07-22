function [V,lam,R] = subspace_iter(V0,A,maxiter,type,a,b,m)

k = size(V0,2);
lam = zeros(k,maxiter);
V = cell(1,maxiter+1);
R = zeros(k,maxiter); 

[V{1},~] = qr(V0,0);
for i = 1:maxiter
    if nargin < 5
        Z = A*V{i};
    else
        Z = chebyshev_filter(A,V{i},m,a,b);
    end
    
    if type == 0
        % orthogonal basis
        [V{i+1},~] = qr(Z,0);
        P = V{i+1}'*A*V{i+1};
        lam(:,i) = diag(P);
        for j = 1:k
            R(j,i) = norm(A*V{i}(:,j) - lam(j,i)*V{i}(:,j)); % residual norm
        end
    else
        % Galerkin problem: generalized eigenvalue of dimension m
        W = Z'*Z;
        P = Z'*A*Z;

        [Y,Lambda] = eig(P,W);
        [lam(:,i),ord] = sort(diag(Lambda),'descend');
        Y = Y(:,ord);

        % Construct Ritz vectors:
        V{i+1} = Z*Y; 
        for j = 1:k
            V{i+1}(:,j) = V{i+1}(:,j)./norm(V{i+1}(:,j));
            R(j,i) = norm(A*V{i}(:,j) - lam(j,i)*V{i}(:,j)); % residual norm
        end
    end
    
    a = max(real(lam(:,i))); % update filter left end-point
    %fprintf('New endpoint: %.2e\n', a)
    
    if mod(i,50) == 0
        fprintf('full iteration: %i \n',i)
    end
end