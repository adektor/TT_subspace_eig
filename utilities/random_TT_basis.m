function V = random_TT_basis(n,L,r,k)

% random TT basis
% n --> dimension at each site
% L --> # of spins
% r --> rank
% k --> dimension of basis
% type --> 'tt' or 'full'

% random basis for subspace
V = cell(1,k);
for i = 1:k
    V{i} = tt_rand(n,L,r); V{i} = round(V{i},1e-12,r);
end