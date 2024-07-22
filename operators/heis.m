function H = heis(L,J)
% Heisenberg hamiltonian
% INPUT: L --> chain length
%        J --> interaction strength
%
% OUTPUT: H --> Hamiltonian (2^L x 2^L)

% pauli matrices
Sx = [0 1; 1 0]; Sy = [0 -1i; 1i 0]; Sz = [1 0; 0 -1];

H = zeros(2^L);
for j = 1:L-1
    H = H + kron3(eye(2^(j-1)),kron(Sx,Sx),eye(2^(L-j-1))) ...
          + kron3(eye(2^(j-1)),kron(Sy,Sy),eye(2^(L-j-1))) ...
          + kron3(eye(2^(j-1)),kron(Sz,Sz),eye(2^(L-j-1))) ...
          ...
          + kron3(eye(2^(j-1)),Sz,eye(2^(L-j)));
end
H = H + kron(eye(2^(L-1)),Sz);

H = -J*H;

end