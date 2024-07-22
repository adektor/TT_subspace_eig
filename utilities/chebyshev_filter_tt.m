function Y = chebyshev_filter_tt(A,X,m,a,b, tol,rmax)

e = (b-a)./2; c = (a+b)./2;
Y = round((round(A*X,tol,rmax) - c*X)/e, tol, rmax);
for i = 2:m
    Yn = round(round(round(A*Y,tol,rmax) - c*Y,tol,rmax)*(2/e) - X, tol, rmax);
    X = Y;
    Y = Yn;
end