function V_lr_full = lr_to_full_basis(V_lr)

% Transform TT basis vectors to full basis vectors. 
% Must ensure TT vectors have small dimension. 

% INPUT: V_lr --> N TT tensors in cell array
% OUTPUT: V_lr_full --> full TT vectors in a M x N matrix, M = n^d

d = V_lr{1}.d;
n = V_lr{1}.n;
maxiter = length(V_lr);

V_lr_full = zeros(n(1)^d,maxiter);

for j = 1:maxiter
    v = full(V_lr{j},V_lr{j}.n');
    V_lr_full(:,j) = reshape(v,[],1);
end

end