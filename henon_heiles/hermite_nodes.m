function [q,L] = hermite_nodes(n)

% calculate Hermite nodes by symmetric companion matrix:
% (Golub/Welsh-type algorithm)
offdiag = sqrt(0.5*[1:n-1]);
comp_matrix = diag(offdiag, -1) + diag(offdiag, 1);
q = sort(eig(comp_matrix));

lambda = 0.11;
%L = -spdiags( [one, -2*one, one], [-1 0 1], n, n) / (h^2);
for i = 1:n
    for j =1:n
        if i == j
            L(i,j) = (4*n - 1 - 2*q(i)^2) / 6;
        else
            L(i,j) = (-1)^(i-j) * (2/(q(i)-q(j))^2 - 0.5);
        end
    end
end

D = spdiags(q, 0, n, n);