function tt = tt_from_cores(C)
% Construct tt_tensor from array of cores C

d = length(C);
r = ones(d+1,1);
n = ones(d,1);
n(1) = size(C{1},3);
core = [];
for k = 1:d
    r(k+1) = size(C{k},2);
    n(k) = size(C{k},3);
    core = [core ; reshape(permute(C{k},[1 3 2]),[],1)] ;
end

tt = tt_tensor();
tt.core = core;
tt.r = r;
tt.d = d;
tt.n = n;
tt.ps = cumsum([1;tt.n.*tt.r(1:tt.d).*tt.r(2:tt.d+1)]) ;