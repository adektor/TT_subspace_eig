function Y = chebyshev_filter_tt(A,X,m,a,b,tol,rmax,ts_mv)

% Chebyshev polynomial filter with low-rank truncations. 
% INPUT: A --> TT matrix
%        X --> TT tensor
%        m --> degree of polyomial
%        a,b --> left,right endpoints for linear map
%        tol --> truncation tolerance
%        rmax --> max. rank
%        var --> 0/1 Compute matvec variationally with AMEn (default 0)

%        OUTPUT: Y --> approximately P(A)X

if tol > 5e-1; tol_mv = 5e-1; end

e = (b-a)./2; c = (a+b)./2;

if ts_mv == 1 % variational matvec
    [AX,~] = axpx(A,X);
else % standard matvec + SVD
    AX = round(A*X,tol,rmax);
end

Y = round((AX - c*X)/e, tol, rmax);
for i = 2:m

    if ts_mv == 1 % variational matvec
        %AY = amen_mv(A,Y,tol_mv,'y0',X);
        AY = axpx(A,Y);
    else % standard matvec + SVD
        AY = round(A*Y,tol,rmax);
    end

    Yn = round(round(AY - c*Y,tol,rmax)*(2/e) - X, tol, rmax);
    X = Y;
    Y = Yn;
end