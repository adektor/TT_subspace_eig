function H = hh_tt(n,d)

% Construct Henon-Heiles TT matrix with mode sizes n and dimension d

% calculate Hermite nodes by symmetric companion matrix:
% (Golub/Welsh-type algorithm)
offdiag = sqrt(0.5*[1:n-1]);
comp_matrix = diag(offdiag, -1) + diag(offdiag, 1);
q = sort(eig(comp_matrix));

%L = -spdiags( [one, -2*one, one], [-1 0 1], n, n) / (h^2);
for i = 1:n
    for j =1:n
        if i == j
            D(i,j) = (4*n - 1 - 2*q(i)^2) / 6;
        else
            D(i,j) = (-1)^(i-j) * (2/(q(i)-q(j))^2 - 0.5);
        end
    end
end
D = reshape(D,[],1); D = permute(D,[3 2 1]);

Q = diag(q); 
Q = reshape(Q,[],1); Q = permute(Q,[3 2 1]);

I = eye(n); I = reshape(I,[],1); I = permute(I,[3 2 1]);

mu = 0.111803;

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
    tt = -0.5*tt_from_cores(cores);
    tt_mat = tt_matrix(tt,n*ones(1,d));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

%% q_k^2
for j = 1:d
    cores = cell(1,d);
    for k = 1:d % cores of rank-1 MPO
        if k == j
            cores{k} = Q.^2;
        else
            cores{k} = I;
        end
    end
    tt = 0.5*tt_from_cores(cores);
    tt_mat = tt_matrix(tt,n*ones(1,d));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

%% q_k^2 q_{k+1}
for j = 1:d-1 % for local operator we make a rank-1 MPO
    cores = cell(1,d);
    for k = 1:d % cores of rank-1 MPO
        if k == j 
            cores{k} = Q.^2;
        elseif k == j+1
            cores{k} = Q;
        else
            cores{k} = I;
        end
    end
    tt = mu*tt_from_cores(cores);
    tt_mat = tt_matrix(tt,n*ones(1,d));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

%% q_{k+1}^3
for j = 1:d-1 % for local operator we make a rank-1 MPO
    cores = cell(1,d);
    for k = 1:d % cores of rank-1 MPO
        if k == j+1
            cores{k} = Q.^3;
        else
            cores{k} = I;
        end
    end
    tt = -(mu./3)*tt_from_cores(cores);
    tt_mat = tt_matrix(tt,n*ones(1,d));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

end