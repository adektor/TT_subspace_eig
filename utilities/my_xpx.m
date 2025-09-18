function y = my_xpx(z,x)

% projects z onto the TT tangent space at x

d = x.d;
xCores = getcores(x);
zCores = getcores(z);

% initialize output cores
yCores = getcores(x);

% orthogonalize y right-to-left and compute right-environments (psi)
psi = cell(d+1,1);
psi{d+1} = 1; psi{1} = 1;

for i = d:-1:2
    [Q,R] = qr_core(xCores{i},'r');
    xCores{i} = Q;
    xCores{i-1} = ncon({xCores{i-1},R},{[-1 1 -3],[1 -2]});
    
    % Compute right environment
    psi{i} = ncon({conj(Q),zCores{i},psi{i+1}}, {[-1 1 2], [-2 3 2], [1 3]});
end

% Now compute the projection while computing left-environments (psi)
for i=1:d

    % Compute cry0 rx(i),rz(i+1),n(i) by contracting 
    % psi{i}   : rx(i),rz(i)
    % zCores{i}: rz(i),rz(i+1),n(i)
    cry0 = ncon({psi{i}, zCores{i}}, {[-1 1],[1 -2 -3]});

    % Compute cry rz(i),rx(i+1),n(i) by contracting
    % cry0     : rz(i),rz(i+1),n(i)
    % psi{i+1} : rx(i+1),rz(i+1)
    cry = ncon({cry0, psi{i+1}}, {[-1 1 -3],[-2 1]});

    % left-orth cry
    [Q,R] = qr_core(cry,'l');
    yCores{i} = Q;
    if i<d
        yCores{i+1} = ncon({R,yCores{i+1}},{[-1 1],[1 -2 -3]});
    end

    % Compute new psi{i+1} rx(i+1),rz(i+1) by contracting 
    % cry0 : rz(i),rz(i+1),n(i)
    % Q    : rz(i),rx(i+1),n(i)
    psi{i+1} = ncon({cry0, conj(Q)}, {[1 -2 2], [1 -1 2]});
end
yCores{d} = R*yCores{d};

y = tt_from_cores(yCores);

end