function V_lr_full = lr_to_full_basis(V_lr)

L = V_lr{1}.d;
maxiter = length(V_lr);

V_lr_full = zeros(2^L,maxiter);

for j = 1:maxiter
    v = full(V_lr{j},V_lr{j}.n');
    V_lr_full(:,j) = reshape(v,[],1);
end

end