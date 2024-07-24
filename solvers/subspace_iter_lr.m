function [V,lam,R,cpu_t,a,b,Y,RQ,gradRQ,PgradRQ] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m,rq,coeff_tol)

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
%        rq --> flag (0 or 1) to compute Rayleigh quotient and grads.
%        coeff_tol --> tolerance for RR coefficients (default 1e-16)

% OUTPUT: V --> TT Ritz vectors (cell array)
%         lam --> approximate (Ritz) values
%         R --> residual of each Ritz vector at each iteration
%         cpu_t --> CPU-time of each iteration
%         a,b, --> updated endpoints for polynomial filter
%         Y --> coefficients of Ritz vectors at each iteration

if nargin == 9 
    coeff_tol = 1e-32;
end

ev_lock=1e-32;

k = length(V0);
lam = zeros(k,maxiter);
V = cell(1,maxiter+1);
R = zeros(k,maxiter); 
Y  = zeros(maxiter,k,k);
RQ = zeros(1,maxiter);
gradRQ = cell(1,maxiter);
PgradRQ = cell(1,maxiter);


V{1} = tt_gs(V0,tol,rmax);
%V{1} = V0;
cpu_t = [];
for i = 1:maxiter
    tic
    Z = cell(1,k);
    for j = 1:k
        if i>2 && abs(R(j,i)-R(j,i-1))>ev_lock
            if nargin < 6
                Z{j} = round(Htt*V{i}{j},tol,rmax);
            else
                Z{j} = chebyshev_filter_tt(Htt,V{i}{j},m,a,b,tol,rmax);
            end
        else
            Z{j} = V{i}{j};
        end
    end

    [V{i+1},lam(:,i),R(:,i),Y(i,:,:)] = RR_lr(Z,Htt,tol,rmax,coeff_tol);
    if rq == 1
        [RQ(i),gradRQ{i},PgradRQ{i}] = rayleigh_quot(V{i}{1},Htt);
    end

    a = max(real(lam(:,i))); % update filter left end-point

    if mod(i,1) == 0
        fprintf('subspace iteration: %i , res. norm: %.2e \n',i,R(1,i))
    end
    cpu_t = [cpu_t, toc];
end