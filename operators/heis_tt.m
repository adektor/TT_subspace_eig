function H = heis_tt(L,J)

% pauli matrices
Sx = [0 1; 1 0];    Sx = reshape(Sx,[],1); Sx = permute(Sx,[3 2 1]);
Sy = [0 -1i; 1i 0]; Sy = reshape(Sy,[],1); Sy = permute(Sy,[3 2 1]);
Sz = [1 0; 0 -1];   Sz = reshape(Sz,[],1); Sz = permute(Sz,[3 2 1]);

I = eye(2); I = reshape(I,[],1); I = permute(I,[3 2 1]);

% initialize Hamiltonian as 0 MPO
ze = tt_zeros(4,L); H = tt_matrix(ze,2*ones(1,L)); 

%% Sxx terms
for j = 1:L-1 % for local operator we make a rank-1 MPO
    cores = cell(1,L);
    for k = 1:L % cores of rank-1 MPO
        if k == j || k == j+1
            cores{k} = Sx;
        else
            cores{k} = I;
        end
    end
    tt = tt_from_cores(cores);
    tt_mat = tt_matrix(tt,2*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

%% Syy terms
for j = 1:L-1 % for local operator we make a rank-1 MPO
    cores = cell(1,L);
    for k = 1:L % cores of rank-1 MPO
        if k == j || k == j+1
            cores{k} = Sy;
        else
            cores{k} = I;
        end
    end
    tt = tt_from_cores(cores);
    tt_mat = tt_matrix(tt,2*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

%% Szz terms
for j = 1:L-1 % for local operator we make a rank-1 MPO
    cores = cell(1,L);
    for k = 1:L % cores of rank-1 MPO
        if k == j || k == j+1
            cores{k} = Sz;
        else
            cores{k} = I;
        end
    end
    tt = tt_from_cores(cores);
    tt_mat = tt_matrix(tt,2*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

%% Sz terms
for j = 1:L % for local operator we make a rank-1 MPO
    cores = cell(1,L);
    for k = 1:L % cores of rank-1 MPO
        if k == j
            cores{k} = Sz;
        else
            cores{k} = I;
        end
    end
    tt = tt_from_cores(cores);
    tt_mat = tt_matrix(tt,2*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

H = -J*H;

end