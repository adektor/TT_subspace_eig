function [eigs2,acc] = process_eigs(eigs,eigs_an)

maxiter = length(eigs);
eigs2 = cell(1,maxiter);
acc = cell(1,maxiter);
for i = 2:maxiter
    eigs2{i} = [];
    for j = 1:i
        eigs2{i} = [eigs2{i} eigs{i}(j)];
        acc{i}(j) = min(abs(eigs{i}(j) - eigs_an));
    end
end

end