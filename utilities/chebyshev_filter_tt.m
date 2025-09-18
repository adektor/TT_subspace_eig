function Y = chebyshev_filter_tt(A,X,k,a,b,tol,rmax,ts)

% Chebyshev polynomial filter with low-rank truncations. 
% INPUT: A --> TT matrix
%        X --> TT tensor
%        m --> degree of polyomial
%        a,b --> left,right endpoints for linear map
%        tol --> truncation tolerance
%        rmax --> max. rank
%        ts --> truncation in tangent space 

%        OUTPUT: Y --> approximately P(A)X

e = (b-a)./2; c = (a+b)./2;

if ts
    I = tt_eye(X.n);
    %Y = ts_proj_matvec((1/e)*(A-c*I), X, X);
    %Y = round(Y, tol, rmax);
    Y = my_axpx((1/e)*(A-c*I), X);
    
    for i = 2:k
        % the following two ts_proj can be combined into 1 using
        % a function that approximates Ax + b on the tangent space at x
        
        %Yn = ts_proj_matvec((2/e)*(A-c*I), Y, Y);
        %Yn = round(Yn, tol, rmax);

        %Yn = axpx((2/e)*(A-c*I), Y);
        
        %Yn = ts_proj_sum({Yn,-1*X},Yn);
        %Yn = round(Yn, tol, rmax);
        
        %Yn = round(Yn-X,tol,rmax);
        
        Yn = my_pax_z(-X,(2/e)*(A-c*I),Y);
        
        X = Y;
        Y = Yn;
    end

else % SVD truncation

    AX = round(A*X,tol,rmax);
    Y = round((AX - c*X)/e, tol, rmax);

    for i = 2:k
        AY = round(A*Y,tol,rmax);
        Yn = round(round(AY - c*Y,tol,rmax)*(2/e) - X, tol, rmax);
        X = Y;
        Y = Yn;
    end

end