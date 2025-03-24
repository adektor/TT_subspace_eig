function gamma = chebyshev_filter_exp(lambda,k,a,b)

% Chebyshev polynomial filter explicit trigonemetric formula. 
% INPUT: lambda --> scalar input (eig. value)
%        k --> degree of polyomial
%        a,b --> left,right endpoints for linear map

%        OUTPUT: gamma --> P(lambda)

e = (b-a)./2; c = (a+b)./2;
lambda = (lambda - c)/e;

if abs(lambda) < 1
    gamma = cos(k*acos(lambda));
else
    gamma = cosh(k*acosh(lambda));
end