function [Q,R] = qr_core(C,s)

% Orthogonalize TT-core C using QR-decomposition
% s = 'l', 'r' for left, right orthogonalization

[r1,r2,n] = size(C);
C = permute(C,[1 3 2]);

if strcmp(s,'r')
    C = reshape(C,[r1,n*r2]);
    C = C.';
    [Q,R] = qr(C,0); 
    R = R.';
    Q = Q.';
    Q = reshape(Q,[r1, n, r2]);
    Q = permute(Q, [1 3 2]);

elseif strcmp(s,'l')
    C = reshape(C,[r1*n,r2]);
    [Q,R]=qr(C,0);
    Q = reshape(Q,[r1, n, r2]);
    Q = permute(Q, [1 3 2]);
end