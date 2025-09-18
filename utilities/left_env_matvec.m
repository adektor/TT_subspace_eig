function L = left_env_matvec(xC,A,vC)

n = A.n; m = A.m;
Att = A.tt;
d = Att.d;
ra = Att.r;
aC = getcores(Att);

L = cell(d-1,1);

cra = reshape(aC{1},[ra(1),ra(2),n(1),m(1)]);
L{1} = ncon({conj(xC{1}),cra,vC{1}}, {[-4 -1 1],[-5 -2 1 2],[-6 -3 2]});

for i = 2:d
    cra = reshape(aC{i},[ra(i),ra(i+1),n(i),m(i)]);
    L{i} = ncon({L{i-1}, conj(xC{i}), cra, vC{i}}, {[1 3 5], [1 -1 2], [3 -2 2 4], [5 -3 4]});
end

end