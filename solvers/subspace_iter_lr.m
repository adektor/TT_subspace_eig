function [V,lam,R,cpu_t,a,b,Y,RQ,gradRQ,PgradRQ,RR_err,T_err] = subspace_iter_lr(V0,Htt,maxiter,tol,rmax,varargin)

% Subspace iteration with low-rank truncations
% INPUT: V0 --> initial subspace of TT tensors (cell array)
%        Htt --> TT matrix
%        maxiter --> number of subspace iterations
%        tol --> truncation tolerance
%        rmax --> max. rank

% OUTPUT: V --> TT Ritz vectors (cell array)
%         lam --> approximate (Ritz) values
%         R --> residual of each Ritz vector at each iteration
%         cpu_t --> CPU-time of each iteration
%         a,b, --> updated endpoints for polynomial filter
%         Y --> coefficients of Ritz vectors at each iteration

%   Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following)
%   The list of option names and default values are:
%        o a - initial left endpoint for polynomial filter 
%        o b - initial right endpoint for polynomial filter 
%        o m - degree of Chebyshev polynomial
%           Omitting a,b,m uses no polymomial filter. Default no filter.
%        o coeff_tol - tolerance for RR coefficients [1e-32]
%        o ev_lock - tolerance for locking eigenvectors [1e-32]
%        o ts_mv - (0 or 1) compute matvecs using tangent space projection [0]

coeff_tol = 1e-100;
ev_lock = 1e-100;

a = [];
b = [];
m = [];
ts_mv = 0;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'a'
            a=varargin{i+1};
        case 'b'
            b=varargin{i+1};
        case 'm'
            m=varargin{i+1};            
        case 'rq'
            rq=varargin{i+1};
        case 'coeff_tol'
            coeff_tol=varargin{i+1};
        case 'ev_lock'
            ev_lock=varargin{i+1};             
        case 'ts_mv'
            ts_mv = varargin{i+1};

        otherwise
            warning('Unrecognized option: %s\n',varargin{i});
    end
end

k = length(V0);
lam = zeros(k,maxiter);
V = cell(1,maxiter+1);
R = zeros(k,maxiter); 
Y  = zeros(maxiter,k,k);
RR_err = zeros(k,maxiter);
RQ = zeros(1,maxiter);
gradRQ = cell(1,maxiter);
PgradRQ = cell(1,maxiter);
Terr = zeros(k,maxiter);

V{1} = tt_gs(V0,tol,rmax);
%V{1} = V0;
cpu_t = [];
for i = 1:maxiter
    tic
    Z = cell(1,k);
    for j = 1:k
        if i>3 && abs(R(j,i-1)-R(j,i-2))>ev_lock
            if isempty(a) % no Chebyshev filter
                if ts_mv == 1 % tangent-space matvec
                    Z{j} = axpx(Htt,V{i}{j});
                else % standard matvec + SVD
                    Z{j} = round(Htt*V{i}{j},tol,rmax);
                end
            else
                Z{j} = chebyshev_filter_tt(Htt,V{i}{j},m,a,b,tol,rmax,ts_mv,R(j,i-1)*1e-2);
                Z{j} = Z{j}./norm(Z{j});
            end
        else
            Z{j} = V{i}{j};
        end
    end

    [V{i+1},lam(:,i),R(:,i),Y(i,:,:),RR_err(:,i)] = RR_lr(Z,Htt,tol,rmax,coeff_tol,ts_mv,tol);

    if mod(i,1) == 0
        cpu_t = [cpu_t, toc];
        Rsort = sort(R(:,i));
        if Rsort(2) < 1e-10
            break;
        end
        fprintf('subspace iteration: %i , res. norm: %.2e \n',i,Rsort(2))
    end
end