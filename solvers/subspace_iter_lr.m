function [V,lam,R,cpu_t,a,b,Y] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m)

% Subspace iteration with low-rank truncations
% INPUT: V0 --> initial subspace of TT tensors (cell array)
%        Htt --> TT matrix
%        maxiter --> number of subspace iterations
%        tol --> truncation tolerance
%        rmax --> max. rank
%        a,b --> initial interval for polynomial filter (updated during
%        iteration) 
%        m --> degree of Chebyshev polynomial
%        Omitting a,b,m uses no polymomial filter.

% OUTPUT: V --> TT Ritz vectors (cell array)
%         lam --> approximate (Ritz) values
%         R --> residual of each Ritz vector at each iteration
%         cpu_t --> CPU-time of each iteration
%         a,b, --> updated endpoints for polynomial filter
%         Y --> coefficients of Ritz vectors at each iteration

k = length(V0);
lam = zeros(k,maxiter);
V = cell(1,maxiter+1);
R = zeros(k,maxiter); 
Y  = zeros(maxiter,k,k);

V{1} = tt_gs(V0,tol,rmax);
%V{1} = V0;
cpu_t = [];
for i = 1:maxiter
    tic
    Z = cell(1,k);
    for j = 1:k
        if nargin < 6
            Z{j} = round(Htt*V{i}{j},tol,rmax);
        else
            Z{j} = chebyshev_filter_tt(Htt,V{i}{j},m,a,b,tol,rmax);
        end
    end

    [V{i+1},lam(:,i),R(:,i),Y(i,:,:)] = RR_lr(Z,Htt,tol,rmax);

    a = max(real(lam(:,i))); % update filter left end-point

    if mod(i,1) == 0
        fprintf('subspace iteration: %i , res. norm: %.2e \n',i,R(1,i))
    end
    cpu_t = [cpu_t, toc];
end