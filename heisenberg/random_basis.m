function V = random_basis(n,L,k)

% random basis
% n --> dimension at each site
% L --> # of spins
% k --> dimension of basis

% random basis for subspace
V = zeros(2^L,k);
for i = 1:k
    V(:,i) = rand(n^L,1);
    V(:,i) = V(:,i)./norm(V(:,i));
end