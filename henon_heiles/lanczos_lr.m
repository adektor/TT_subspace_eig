function [V,T,eigs] = lanczos_lr(Att,v0tt,maxiter,tol,rmax,a,b,m)

% Approximate Lanczos with TT vectors/matrix 
% 
% INPUT: Att --> tt-matrix
%        v0tt --> tt-vector (norm 1)
%        maxiter --> number of Lanczos steps
%        rmax --> max tt-rank
%        tol --> truncation tolerance
% OUTPUT: V --> cell array containing tt-vectors

alpha = [];
beta = 0;
eigs = cell(1,maxiter);

V{1} = v0tt;

if nargin < 6
    V{2} = round(Att*V{1},tol,rmax);
else
    V{2} = chebyshev_filter_tt(Att,V{1},m,a,b,tol,rmax);
end

alpha(1) = dot(V{2},V{1});

% two-term recurrence
V{2} = round( V{2} - alpha(1)*V{1},tol,rmax);

W = [dot(V{1},V{1})];
P = dot(V{1},Att*V{1});

for j = 2:maxiter
    beta = [beta norm(V{j})];
    V{j} = V{j}./beta(j);

    % new Krylov vector
    if nargin < 6
        V{j+1} = round(Att*V{j},tol,rmax);
    else
        V{j+1} = chebyshev_filter_tt(Att,V{j},m,a,b,tol,rmax);
    end
    %V{j+1} = round(Att*V{j},tol,rmax);

    alpha = [alpha dot(V{j+1},V{j})];

    % three-term recurrence
    V{j+1} = round(V{j+1} - alpha(j)*V{j},  tol, rmax);
    V{j+1} = round(V{j+1} - beta(j)*V{j-1}, tol, rmax);

    % add new row and column to overlap matrices: 
    for k = 1:j
        W(j,k) = dot(V{j},V{k});
        W(k,j) = W(j,k);

        P(j,k) = dot(V{j},Att*V{k});
        P(k,j) = P(j,k);
    end

    eigs{j} = eig(P,W); 
    eigs{j} = sort(real(eigs{j}),'descend');
    
    if mod(j,10) == 0
        fprintf("Lanczos iteration %i: \t Condition # of overlap: %.3e \n",j, cond(W))
    end
end

T = diag(alpha) + diag(beta(2:end),-1) + diag(beta(2:end),1);

end