function [RQ,gradRQ,PgradRQ] = rayleigh_quot(x,A)

% Compute the Rayleigh-Quotient of x with matrix A
% INPUTS: 
% x --> TT tensor or full vector
% A --> TT matrix or dense matrix 

% OUTPUTS: 
% RQ --> Rayleigh quotient
% gradRQ --> gradient of RQ
% PgradRQ --> projection of gradient onto tangent space

x = x./norm(x);
RQ = dot(x,A*x);

gradRQ = 2*(A*x - RQ*x);

if isa(x,'tt_tensor')
    PgradRQ = 2*(axpx(A,x) - RQ*x);
else
    PgradRQ = 0;
end

end