function y = my_pax_z(z,a,x)

% projects Ax + z onto the TT tangent space at x

d = x.d;
xCores = getcores(x);
zCores = getcores(z);

n = a.n; m = a.m;
att = a.tt;
ra = att.r;
aCores = getcores(att);

% initialize output cores
yCores = getcores(x);

% orthogonalize y right-to-left and compute right-environments (psi)
psi = cell(d+1,2);
psi{d+1,1} = 1; psi{1,1} = 1;
psi{d+1,2} = 1; psi{1,2} = 1;

for i = d:-1:2
    [Q,R] = qr_core(xCores{i},'r');
    xCores{i} = Q;
    xCores{i-1} = ncon({xCores{i-1},R},{[-1 1 -3],[1 -2]});
    
    % right environment for z
    psi{i,1} = ncon({conj(Q),zCores{i},psi{i+1,1}}, {[-1 1 2], [-2 3 2], [1 3]});

    % right environment for Av
    cra = reshape(aCores{i},[ra(i),ra(i+1),n(i),m(i)]);

    % Compute psi{i} ra(i),rx(i),ry(i) by contracting
    % psi{i+1}  : ra(i+1),rx(i+1),ry(i+1)
    % xCores{i} : rx(i),rx(i+1),m(i)
    % cra       : ra(i),ra(i+1),n(i),m(i)
    % yCores{i} : ry(i),ry(i+1),n(i)
    psi{i,2} = ncon({xCores{i},psi{i+1,2}}, {[-1 1 -2],[-3 1 -4]}); % rx(i),m(i),ra(i+1),ry(i+1)
    psi{i,2} = ncon({cra,psi{i,2}}, {[-3 1 -4 2],[-1 2 1 -2]});     % rx(i),ry(i+1),ra(i),n(i)
    psi{i,2} = ncon({conj(Q),psi{i,2}}, {[-3 1 2],[-2 1 -1 2]});    % ra(i),rx(i),ry(i)
end

% Now compute the projection while computing left-environments (psi)
for i=1:d
    cra = reshape(aCores{i},[ra(i),ra(i+1),n(i),m(i)]);
    cry0A = ncon({psi{i,2}, cra, xCores{i}}, {[2 3 -1],[2 -3 -2 1],[3 -4 1]});
    cryA = ncon({cry0A, psi{i+1,2}}, {[-1 -3 1 2],[1 2 -2]});

    % Compute cry0 rx(i),rz(i+1),n(i) by contracting 
    % psi{i}   : rx(i),rz(i)
    % zCores{i}: rz(i),rz(i+1),n(i)
    cry0 = ncon({psi{i,1}, zCores{i}}, {[-1 1],[1 -2 -3]});

    % Compute cry rz(i),rx(i+1),n(i) by contracting
    % cry0     : rz(i),rz(i+1),n(i)
    % psi{i+1} : rx(i+1),rz(i+1)
    cry = ncon({cry0, psi{i+1,1}}, {[-1 1 -3],[-2 1]});

    % left-orth cry
    [Q,R] = qr_core(cry+cryA,'l');
    yCores{i} = Q;
    if i<d
        yCores{i+1} = ncon({R,yCores{i+1}},{[-1 1],[1 -2 -3]});
    end

    % left-environement psi{i+1,1} rx(i+1),rz(i+1) for z by contracting 
    % cry0 : rz(i),rz(i+1),n(i)
    % Q    : rz(i),rx(i+1),n(i)
    psi{i+1,1} = ncon({cry0, conj(Q)}, {[1 -2 2], [1 -1 2]});

    % left-environement psi{i+1,2} for Av
    psi{i+1,2} = ncon({cry0A, conj(Q)},{[2 1 -1 -2], [2 -3 1]});
end
yCores{d} = R*yCores{d};

y = tt_from_cores(yCores);

end