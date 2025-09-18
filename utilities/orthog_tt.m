function C = orthog_tt(C,r,s)
% Orthogonalize a TT either from left or right. 
%
% INPUT:
% C --> TT cores
% r --> TT rank
% s --> string containing 'l' or 'r'
%
% OUTPUT: 
% C --> orthogonalized TT cores

d = length(C);

if s == 'l'
    for i = 1:d-1
        n = size(C{i},3); % number of points in core i
        A = permute(C{i},[1 3 2]); % move physical dimension to middle
        A = reshape(A,r(i)*n,r(i+1));
        [U,S,V] = svd(A,'econ');
        r(i+1) = size(U,2); % new rank
        C{i} = reshape(U,r(i),n,r(i+1)); % core i is now orthogonal
        C{i} = permute(C{i},[1 3 2]); 
        C{i} = reshape(C{i},r(i),r(i+1),n);

        % move S*V' into core i+1
        C{i+1} = ncon({(S*V'),C{i+1}}, {[-1,1],[1,-2,-3]});
    end
    
elseif s == 'r'
    for i = d:-1:2
        n = size(C{i},3); % number of points in core i
        A = permute(C{i},[1 3 2]); % move physical dimension to middle
        A = reshape(A,r(i),n*r(i+1));
        [U,S,V] = svd(A,'econ');
        r(i) = size(U,2); % new rank
        C{i} = reshape(V',r(i),n,r(i+1)); % core i is now orthogonal
        C{i} = permute(C{i},[1 3 2]); 
        C{i} = reshape(C{i},r(i),r(i+1),n);

        % move S*U into core i-1
        C{i-1} = ncon({C{i-1},(U*S)}, {[-1,1,-3],[1,-2]});
    end
end

end