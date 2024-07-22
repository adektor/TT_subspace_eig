function [lam,psi] = exact_eigs(Htt)
    
    H = full(Htt);                          
    [psi,lam] = eig(H);
    [lam,ord] = sort(real(diag(lam)),'ascend');
    psi = psi(:,ord);

end