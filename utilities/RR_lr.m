function [V,lam,r,Y,RR_err] = RR_lr(Z,Htt,tol,rmax,coeff_tol,var_mv,tol_mv)

% Rayleigh-Ritz projection with non-orthogonal low-rank basis vectors 
% by solving a generalized eigenvalue problem
%
% INPUT: Z --> low-rank basis vectors stored in cell array
%        Htt --> MPO
%        tol --> truncation tolerance
%        rmax --> max. rank
%        coeff_tol --> tolerance for RR coefficients 
%
% OUTPUT: V --> low-rank Ritz vectors
%         lam --> Ritz values
%         R --> residual norms
%         Y --> coefficient vectors from generalized eigenvalue problem

k = length(Z); % dimension of subspace

AZ = cell(1,k);
RR_err = zeros(k,1);
for j = 1:k
    if var_mv == 1 % variational matvec
        [AZ{j},~] = amen_mv(Htt,Z{j},tol_mv,'y0',Z{j});
    else % standard matvec + SVD
        AZ{j} = round(Htt*Z{j},tol,rmax);
        RR_err(j) = norm(AZ{j} - Htt*Z{j})./norm(Htt*Z{j});
    end
end

% Galerkin problem: generalized eigenvalue of dimension m
W = overlap_mat(Z,Z);
P = overlap_mat(Z,AZ);

[Y,Lambda] = eig(P,W);
Lambda = real(Lambda);
[lam,ord] = sort(diag(Lambda),'descend');
Y = Y(:,ord);

V = cell(1,k); r = zeros(1,k);

% Construct Ritz vectors:
for j = 1:k
    %V{j} = Y(1,j)*Z{1};
    [~,ind] = max(abs(Y(:,j)));
    V{j} = Y(ind,j)*Z{ind};
    
    for p = 2:k
        if p ~= ind && abs(Y(p,j))>coeff_tol
            V{j} = round(V{j} + Y(p,j)*Z{p},tol,rmax);
        end
    end
    V{j} = V{j}./norm(V{j});

    if var_mv == 1 % variational matvec
        r(j) = norm(amen_mv(Htt,V{j},tol_mv,'y0',V{j}) - lam(j)*V{j}); % residual norm
    else % standard matvec + SVD
        r(j) = norm(Htt*V{j} - lam(j)*V{j}); % residual norm
    end
end

end