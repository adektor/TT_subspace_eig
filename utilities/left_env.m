function L = left_env(C,dC)
d = length(C);
L = cell(d-1,1);
L{1} = ncon({conj(C{1}),dC{1}}, {[-3 -1 1],[-4 -2 1]});
for i = 2:d-1
    %L{i} = ncon({L{i-1},dC{i}}, {[-1 1],[1 -2 -3]});
    %L{i} = ncon({L{i},conj(C{i})}, {[1 -2 2],[1 -1 2]});
    L{i} = ncon({L{i-1}, conj(C{i}), dC{i}}, {[1 3],[1 -1 2],[3 -2 2]});
end

end