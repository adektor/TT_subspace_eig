function H = heis_tt_spin1(L)

% TT matrix for Spin-1 Heisenberg chain of length L

% pauli matrices
Sx = (1/sqrt(2))*[0 1 0; 1 0 1; 0 1 0];         Sx = reshape(Sx,[],1); Sx = permute(Sx,[3 2 1]);
Sy = (1/sqrt(2))*[0 -1i 0; 1i 0 -1i; 0 1i 0];   Sy = reshape(Sy,[],1); Sy = permute(Sy,[3 2 1]);
Sz = [1 0 0; 0 0 0; 0 0 -1];                    Sz = reshape(Sz,[],1); Sz = permute(Sz,[3 2 1]);

I = eye(3); I = reshape(I,[],1);                I = permute(I,[3 2 1]);

% initialize Hamiltonian as 0 MPO
ze = tt_zeros(9,L); H = tt_matrix(ze,3*ones(1,L)); 

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
    tt_mat = tt_matrix(tt,3*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

% 'wrap around' (periodic) term
cores = cell(1,L);
for k = 1:L
    cores{k} = I;
end
cores{L} = Sx; 
cores{1} = Sx;
tt = tt_from_cores(cores);
tt_mat = tt_matrix(tt,3*ones(1,L));
H = H + tt_mat; % add rank-1 MPO to Hamiltonian

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
    tt_mat = tt_matrix(tt,3*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

% 'wrap around' (periodic) term
cores = cell(1,L);
for k = 1:L
    cores{k} = I;
end
cores{L} = Sy; 
cores{1} = Sy;
tt = tt_from_cores(cores);
tt_mat = tt_matrix(tt,3*ones(1,L));
H = H + tt_mat; % add rank-1 MPO to Hamiltonian

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
    tt_mat = tt_matrix(tt,3*ones(1,L));
    H = H + tt_mat; % add rank-1 MPO to Hamiltonian
end

% 'wrap around' (periodic) term
cores = cell(1,L);
for k = 1:L
    cores{k} = I;
end
cores{L} = Sz; 
cores{1} = Sz;
tt = tt_from_cores(cores);
tt_mat = tt_matrix(tt,3*ones(1,L));
H = H + tt_mat; % add rank-1 MPO to Hamiltonian

end