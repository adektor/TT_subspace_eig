function [lam,psi] = exact_eigs(Htt)
% Transform a TT matrix to full matrix and run eigs. 
    
    H = full(Htt);                          
    [psi,lam] = eig(H);
    [lam,ord] = sort(real(diag(lam)),'ascend');
    psi = psi(:,ord);

end