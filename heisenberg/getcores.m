function G = getcores(tt)
% Input: tt --> tt_tensor
% Output: G --> cell array containing the d tensor cores of tt
%               
% Each cell of G contains a 3 tensor with dim. ordering
% 
%    1   [-------]   2
%    ----|  G{i} |----
%        [-------]
%            |
%            | 3


G = cell(1,tt.d);

d = tt.d;
core = tt.core;
ps = tt.ps;
r = tt.r;
M = tt.n;

for j = 1:d
    G{j} = core(ps(j):ps(j+1)-1);
    G{j} = reshape(G{j},r(j),M(j),r(j+1));
    G{j} = permute(G{j},[1 3 2]);
end