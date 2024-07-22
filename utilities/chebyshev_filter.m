function Y = chebyshev_filter(A,X,m,a,b)

% Chebyshev polynomial filter. 
% INPUT: A --> matrix
%        X --> vector
%        m --> degree of polyomial
%        a,b --> left,right endpoints for linear map

%        OUTPUT: Y --> P(A)X


e = (b-a)./2; c = (a+b)./2;
Y = (A*X - c*X)/e;
for i = 2:m
    Yn = (A*Y - c*Y)*(2/e) - X;
    X = Y;
    Y = Yn;
end