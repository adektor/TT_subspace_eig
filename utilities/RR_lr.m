function [V,lam,r] = RR_lr(Z,A,tol,rmax,ts)

% Rayleigh-Ritz projection with non-orthogonal low-rank basis vectors 
% by solving a generalized eigenvalue problem
%
% INPUT: Z --> low-rank basis vectors stored in cell array
%        Htt --> MPO
%        tol --> truncation tolerance
%        rmax --> max. rank
%        coeff_tol --> tolerance for RR coefficients 
%        ts --> truncation in tangent space 
%
% OUTPUT: V --> low-rank Ritz vectors
%         lam --> Ritz values
%         R --> residual norms
%         Y --> coefficient vectors from generalized eigenvalue problem

m = length(Z); % dimension of subspace

AZ = cell(1,m);
for j = 1:m
    if ts % tangent-space matvec 
        AZ{j} = axpx(A,Z{j});
    else % matvec + SVD
        AZ{j} = round(A*Z{j},tol,rmax);
    end
end

% Generalized m x m eigenvalue problem
W = overlap_mat(Z,Z);
P = overlap_mat(Z,AZ);

[Y,Lambda] = eig(P,W);
Lambda = real(Lambda);
[lam,ord] = sort(diag(Lambda),'descend');
Y = Y(:,ord);

V = cell(1,m); r = zeros(1,m);

% Construct Ritz vectors
for j = 1:m
    % V_j = \sum_{p=1}^m Y(p,j)*Z_p
    if ts % tangent space truncation
        YZ = cellfun(@(tensor, scalar) scalar * tensor, Z', num2cell(Y(:,j)), 'UniformOutput', false);
        V{j} = ts_proj_sum(YZ,Z{j});
        V{j} = round(V{j},tol,rmax);

    else % SVD
        V{j} = Y(1,j)*Z{1};
        [~,ind] = max(abs(Y(:,j)));
        V{j} = Y(ind,j)*Z{ind};
    
        for p = 2:m
            if p ~= ind
                V{j} = round(V{j} + Y(p,j)*Z{p},tol,rmax);
            end
        end
    end

    V{j} = V{j}./norm(V{j});
    r(j) = norm(axpx(A,V{j}) - lam(j)*V{j}); % residual norm (cpu intensive)
end

end