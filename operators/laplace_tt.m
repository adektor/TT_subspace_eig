function H = laplace_tt(n,d)

% finite difference stencil
one = ones(n,1);
D = -2*diag(one) + diag(one(1:end-1),-1) + diag(one(1:end-1),1);
D = reshape(D,[],1); D = permute(D,[3 2 1]);

I = eye(n); I = reshape(I,[],1); I = permute(I,[3 2 1]);

% initialize Hamiltonian as 0 MPO
ze = tt_zeros(n^2,d); H = tt_matrix(ze,n*ones(1,d));

%% Laplace
for j = 1:d
    cores = cell(1,d);
    for k = 1:d % cores of rank-1 MPO
        if k == j
            cores{k} = D;
        else
            cores{k} = I;
        end
    end
    tt = tt_from_cores(cores);
    tt_mat = tt_matrix(tt,n*ones(1,d));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

end