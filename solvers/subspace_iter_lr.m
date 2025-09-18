function [V,lam,R,cpu_t,a,b] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,a,b,m,ts,lam_max,verb)

% Subspace iteration with low-rank truncations
% INPUT: V0 --> initial subspace of TT tensors (cell array)
%        Htt --> TT matrix
%        maxiter --> number of subspace iterations
%        tol --> truncation tolerance
%        rmax --> max. rank
%        a --> initial left endpoint for polynomial filter 
%        b --> initial right endpoint for polynomial filter 
%        m --> degree of Chebyshev polynomial
%        ts --> 0 or 1, perform rank-truncation in tangent space
%        lam_max --> max eig estimate used for shift
%        verb --> 0 silent, 1 print leading eigenvalues

% OUTPUT: V --> TT Ritz vectors (cell array)
%         lam --> approximate (Ritz) values
%         R --> residual of each Ritz vector at each iteration
%         cpu_t --> CPU-time of each iteration
%         a,b, --> updated endpoints for polynomial filter


if nargin < 10 % no shift
    lam_max = 0;
end

coeff_tol = 1e-100;

k = length(V0);
lam = zeros(k,maxiter);
V = cell(1,maxiter+1);
R = zeros(k,maxiter); 

V{1} = tt_gs(V0,tol,rmax);
%V{1} = V0;
cpu_t = [];
for i = 1:maxiter
    tic
    Z = cell(1,k);
    for j = 1:k
        if isempty(a) % no Chebyshev filter
            if ts == 1 % tangent-space matvec
                Z{j} = axpx(Htt,V{i}{j});
            else % standard matvec + SVD
                Z{j} = round(Htt*V{i}{j},tol,rmax);
            end
        else
            Z{j} = chebyshev_filter_tt(Htt,V{i}{j},m,a,b,tol,rmax,ts);
            Z{j} = Z{j}./norm(Z{j});
        end
    end

    [V{i+1},lam(:,i),R(:,i)] = RR_lr(Z,Htt,tol,rmax,coeff_tol);

    if verb > 0
        cpu_t = [cpu_t, toc];
        % Rsort = sort(R(:,i));
        % if Rsort(2) < 1e-10
        %     break;
        % end
        esort = sort(lam(:,i),'ascend');
        fprintf('subspace iteration: %i | energy 1: %.8f | energy 2: %.8f | runtime %.f min \n',i,esort(1)+lam_max,esort(2)+lam_max, cpu_t(end)/60);
        % fprintf('subspace iteration: %i , res. norm: %.2e \n',i,Rsort(1))
    end
end