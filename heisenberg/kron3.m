function K = kron3(A,B,C)
%KRON3  Kronecker tensor product.

K = kron(A,kron(B,C));

end
