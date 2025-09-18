function R = right_env_matvec(xC,A,vC)
d = length(xC);

n = A.n; m = A.m;
Att = A.tt;
ra = Att.r;
aC = getcores(Att);

R = cell(d-1,1);

cra = reshape(aC{d},[ra(d),ra(d+1),n(d),m(d)]);
R{d} = ncon({conj(xC{d}),cra,vC{d}}, {[-1 -4 1],[-2 -5 1 2],[-3 -6 2]});

for i = d-1:-1:2
    cra = reshape(aC{i},[ra(i),ra(i+1),n(i),m(i)]);
    R{i} = ncon({conj(xC{i}),cra,vC{i},R{i+1}}, {[-1 1 2],[-2 3 2 4],[-3 5 4],[1 3 5]});
end

end