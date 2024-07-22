% Adapted from:
%   TTeMPS Toolbox. 
%   Michael Steinlechner, 2013-2016
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

function crx = shift_block_core( crx, mu, nu, tol, maxrank )
% INPUTS: 
%   crx --> TT cores
%   mu --> current block core
%   nu --> nu=mu +/- 1 desired block core (left or right of mu)
%   p --> # of blocks

r1 = size(crx{mu},1);
n = size(crx{mu},2);
r2 = size(crx{mu},3);
p = size(crx{mu},4);

    if ~exist('tol', 'var')
        tol = eps;
    end
    if ~exist('maxrank', 'var')
        maxrank = inf;%r(mu+1);
    end

    if mu == nu-1   
        % shift block one to the right
        U = permute( crx{mu}, [1, 2, 4, 3] );
        U = reshape( U, [r(mu)*n(mu), p*r(mu+1)] );

        [U,S,V] = svd( U, 'econ' );
        if p == 1 
            s = length(diag(S));
        else
            s = trunc_singular( diag(S), tol );
        end
        if length(diag(S)) >= s+1
            disp(['cut singular value of rel. magnitude (s_{i+1}/s_1): ', ...
                        num2str(S(s+1,s+1)/S(1,1))])
        end
        U = U(:,1:s);
        crx{mu} = reshape( U, [r(mu), n(mu), s] );
        W = S(1:s,1:s)*V(:,1:s)';
        W = reshape( W, [s, p, r(mu+1)]);
        W = permute( W, [1, 3, 2]);
        
        C = zeros( [s, n(nu), r(nu+1), p] ); 
        for k = 1:p
            C(:,:,:,k) = tensorprod_ttemps( crx{nu}, W(:,:,k), 1);
        end
        
        crx{nu} = C;

    elseif mu == nu+1
        % shift block one to the left
        V = permute( crx{mu}, [1, 4, 2, 3] );
        V = reshape( V, [r1*p, n*r2] );

        [U,S,V] = svd( V, 'econ' );
        s = my_chop2(diag(S),tol);

        if length(diag(S)) >= s+1
            disp(['cut singular value of rel. magnitude (s_{i+1}/s_1): ', ...
                        num2str(S(s+1,s+1)/S(1,1))])
        end

        V = V(:,1:s)';
        crx{mu} = reshape( V, [s, n, r2] );

        W = U(:,1:s)*S(1:s,1:s);
        W = reshape( W, [r1, p, s]);
        W = permute( W, [1, 3, 2]);

        C = zeros([size(crx{nu},1,2),s,p]);
        for k = 1:p
            C(:,:,:,k) = ncon({crx{nu},W(:,:,k)},{[-1 -2 1],[1 -3]});
        end
        
        crx{nu} = C;
    else
        error('Can only shift the superblock one core left or right')
    end

end