function Q = tt_gs(V,tol,rmax)

% Approximate Gram-Schmidt for TT tensors. 
% SVD rounding is applied after each sum. 

% INPUT:    V --> basis of TT_tensors (cell array)
%           tol --> truncation tolerance
%           rmax --? max. rank

% OUTPUT:   Q --> Approximately orthonormalized basis of TT tensors

k = length(V);
Q = cell(1,k);

Q{1} = V{1}./norm(V{1});

for i = 2:k
    Q{i} = V{i};
    for j = 1:i-1
        Q{i} = round(Q{i} - dot(conj(Q{i}),Q{j})*Q{j} ,tol,rmax);
        %dt = dot(Q{i},Q{j})
    end
    Q{i} = Q{i}./norm(Q{i});
end

end