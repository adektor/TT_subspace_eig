function [W] = overlap_mat(V1,V2)

% Construct overlap matrix
% INPUT: V --> cell array containing basis of N TT vectors
% OUTPUT: N x N matrix W_ij = < V_i | V_j >


N = length(V1); % dimension of basis
W = zeros(N,N);
for i = 1:N
    for j = 1:N
        W(i,j) = dot(V1{i},V2{j});
    end
end

end