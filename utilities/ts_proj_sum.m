function Pz = ts_proj_sum(z,x)
% Project \sum_p z{p} onto the TT tangent space at x
% z{p} and x are tt_tensors
%
% First implementation (not optimized)
% compute all left/right orthogonalized cores first
% ALSO JUST COMPUTE ALL RIGHT/LEFT ENVS ONCE IN THE BEGINNING!

d = x.d;
xC = getcores(x);
xCl = orthog_tt(xC,x.r,'l');
xCr = orthog_tt(xC,x.r,'r');

dC = cell(d,1);  % stores variations of cores
PzC = cell(d,1); % stores TT-cores of Pz

Ns = length(z); % number of summands
zC = cell(d,1); % stores cores of each summand

L = cell(Ns,1); R = cell(Ns,1);
for p = 1:Ns
    zC{p} = getcores(z{p});
    L{p} = left_env(xCl,zC{p});
    R{p} = right_env(xCr,zC{p});
end

% first core
dC{1} = zeros(size(xC{1}));
for p = 1:Ns
    tmp = ncon({zC{p}{1},R{p}{2}}, {[-1 2 -3],[-2 2]});
    tmp = tmp - ncon({xCl{1},conj(xCl{1}),tmp}, {[-1 3 -3],[2 3 1],[2 -2 1]});
    dC{1} = dC{1} + tmp;
end
PzC{1} = cat(2,dC{1},xCl{1});

for i = 2:d-1 % interior cores
    dC{i} = zeros(size(xC{i}));
    for p = 1:Ns
        tmp = ncon({L{p}{i-1},zC{p}{i},R{p}{i+1}}, {[-1 1],[1 2 -3],[-2 2]});
        tmp = tmp - ncon({xCl{i},conj(xCl{i}),tmp}, {[-1 3 -3],[2 3 1],[2 -2 1]});
        dC{i} = dC{i} + tmp;
    end
    zeroblock = zeros(size(xCr{i}));
    tmp1 = cat(2,xCr{i},zeroblock);
    tmp2 = cat(2,dC{i},xCl{i});
    PzC{i} = cat(1,tmp1,tmp2);
end

% final core
dC{d} = zeros(size(xC{d}));
for p = 1:Ns
    dC{d} = dC{d} + ncon({L{p}{d-1},zC{p}{d}}, {[-1 1],[1 -2 -3]});
end
PzC{d} = cat(1,xCr{d},dC{d});

Pz = tt_from_cores(PzC);

end