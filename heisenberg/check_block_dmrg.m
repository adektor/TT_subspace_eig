function [R_dmrg,numel_blk] = check_block_dmrg(x,theta,Htt,tol)

% Move block around, check residual norm and memory
%
% INPUT: x --> tt tensor w/ block in last core
%        theta --> approximate eigenvalues from block DMRG
%        Htt --> TT-matrix
%        tol --> truncation tolerance for moving block core around
% 
% OUTPUT: R_dmrg --> residual norm. Row index specifies which core has
%                    block, column index specifies eigenvector 
%         numel_blk --> total # of elements in the block TT tensor 
%                       for different block cores. 


% Extract eigenvectors from block-TT, compute residual norm
d = x.d;
r = x.r;
n = x.n;

crx = core2cell(x);
k = size(crx{end},3);
crx{end} = reshape(crx{end},r(end-1),n(d),1,k);

tt = cell(1,k);
p=1;
R_dmrg = zeros(d-1,k);
numel_blk = zeros(1,d);

for i = 1:k
    Ci = crx(1:d-1);
    Ci{end+1} = crx{d}(:,:,:,i);
    tt{i} = cell2core(tt_tensor,Ci);
    tt{i} = round(tt{i},1e-12,1024);
    R_dmrg(p,i) = norm(Htt*tt{i}-theta(i)*tt{i});
end

% count elements
for i = 1:d
    numel_blk(p) = numel_blk(p) + numel(crx{i});
end

for sft = d:-1:2
    p = p + 1;
    crx = shift_block_core( crx, sft, sft-1,tol);
    for i = 1:k
        Ci = cell(1,d);
        Ci(1:sft-2) = crx(1:sft-2);
        Ci{sft-1} = crx{sft-1}(:,:,:,i);
        Ci(sft:end) = crx(sft:end);
        tt{i} = cell2core(tt_tensor,Ci);
        tt{i} = round(tt{i},1e-12,1024);
        R_dmrg(p,i) = norm(Htt*tt{i}-theta(i)*tt{i});
    end
    for i = 1:d
    numel_blk(p) = numel_blk(p) + numel(crx{i}); 
    end
end
