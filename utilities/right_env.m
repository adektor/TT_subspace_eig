function R = right_env(C,dC)
d = length(C);
R = cell(d-1,1);
R{d} = ncon({conj(C{end}),dC{end}}, {[-1 -3 1],[-2 -4 1]});
for i = d-1:-1:2
    %R{i} = ncon({R{i+1},dC{i}}, {[-1 1],[-2 1 -3]});
    %R{i} = ncon({R{i},conj(C{i})}, {[1 -2 2],[-1 1 2]});
    R{i} = ncon({conj(C{i}), dC{i}, R{i+1}}, {[-1 1 2], [-2 3 2], [1 3]});
end

end