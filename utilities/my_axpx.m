function [y,psi] = my_axpx(a,x)

% TT-matrix vector product Ax and then projects the 
% result onto TT tangent space at x

d = x.d;
xCores = getcores(x);

n = a.n; m = a.m;
att = a.tt;
ra = att.r;
aCores = getcores(att);

% initialize output
y = x;
yCores = getcores(y);

% orthogonalize y right-to-left and compute right-environments (psi)
psi = cell(d+1,1);
psi{d+1} = 1; psi{1} = 1;

for i = d:-1:2
    [Q,R] = qr_core(yCores{i},'r');
    yCores{i} = Q;
    yCores{i-1} = ncon({yCores{i-1},R},{[-1 1 -3],[1 -2]});
    
    cra = reshape(aCores{i},[ra(i),ra(i+1),n(i),m(i)]);
    % Compute psi{i} ra(i),rx(i),ry(i) by contracting
    % psi{i+1}  : ra(i+1),rx(i+1),ry(i+1)
    % xCores{i} : rx(i),rx(i+1),m(i)
    % cra       : ra(i),ra(i+1),n(i),m(i)
    % yCores{i} : ry(i),ry(i+1),n(i)
    psi{i} = ncon({xCores{i},psi{i+1}}, {[-1 1 -2],[-3 1 -4]}); % rx(i),m(i),ra(i+1),ry(i+1)
    psi{i} = ncon({cra,psi{i}}, {[-3 1 -4 2],[-1 2 1 -2]});     % rx(i),ry(i+1),ra(i),n(i)
    psi{i} = ncon({conj(Q),psi{i}}, {[-3 1 2],[-2 1 -1 2]});    % ra(i),rx(i),ry(i)
end

% Now compute the projection while computing left-environments (psi)
for i=1:d
    cra = reshape(aCores{i},[ra(i),ra(i+1),n(i),m(i)]);
    % Compute cry0 ry(i),n(i),ra(i+1),rx(i+1) by contracting 
    % psi{i}   : ra(i),rx(i),ry(i)
    % cra      : ra(i),ra(i+1),n(i),m(i)
    % xCores{i}: rx(i),rx(i+1),m(i)
    cry0 = ncon({psi{i}, cra, xCores{i}}, {[2 3 -1],[2 -3 -2 1],[3 -4 1]});

    % Compute cry ry(i),ry(i+1),n(i) by contracting
    % cry0     : ry(i),n(i),ra(i+1),rx(i+1)
    % psi{i+1} : ra(i+1),rx(i+1),ry(i+1)
    cry = ncon({cry0, psi{i+1}}, {[-1 -3 1 2],[1 2 -2]});

    % left-orth cry
    [Q,R] = qr_core(cry,'l');
    yCores{i} = Q;
    if i<d
        yCores{i+1} = ncon({R,yCores{i+1}},{[-1 1],[1 -2 -3]});
    end

    % Compute new psi{i+1} ra(i+1),rx(i+1),ry(i+1) by contracting 
    % cry0 : ry(i),n(i),ra(i+1),rx(i+1)
    % Q    : ry(i),ry(i+1),n(i)
    psi{i+1} = ncon({cry0, conj(Q)},{[2 1 -1 -2], [2 -3 1]});
end
yCores{d} = R*yCores{d};

y = tt_from_cores(yCores);

end