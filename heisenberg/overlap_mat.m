function [W] = overlap_mat(V1,V2)
% Construct reorthogonalization matrix
% INPUT: V --> cell array containing basis of TT vectors

N = length(V1); % dimension of basis

% Construct overlap matrix W_ij = < V_i | V_j >
W = zeros(N,N);
for i = 1:N
    for j = 1:N
        W(i,j) = dot(V1{i},V2{j});
    end
end

end