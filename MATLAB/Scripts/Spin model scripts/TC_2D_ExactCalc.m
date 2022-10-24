Dim = [3 3]; N = prod(Dim); Nv = 2*N; ID = speye(2^Nv);

GraphObj = SquareDualGraph(HypCub(Dim,[1 1],eye(2),0));
VertList = GraphObj.ExtraLabels{1}; PlaqList = GraphObj.ExtraLabels{2}; z = 4;

basis = zeros(2^Nv,Nv);
for k=1:2^Nv
  v = ind2dig(k,2*ones(1,Nv));
  basis(k,:) = 3-2*v; % Express basis as sz = {+1,-1}
end

[Sx,~,Sz,~,~] = spin_half_operators(Nv);

% Vertex contributions.
SVT = sparse(1,1,0,2^Nv,2^Nv,2^Nv);
for v = 1:N
    SV = speye(2^Nv);
    for s = 1:z
        SV = SV * Sz{VertList(v,s)};
    end
    SVT = SVT + SV;
end
% Plaquette contributions.
SPT = sparse(1,1,0,2^Nv,2^Nv,2^Nv);
for p = 1:N
    SP = speye(2^Nv);
    for s = 1:z
        SP = SP * Sx{PlaqList(p,s)};
    end
    SPT = SPT + SP;
end

H = - SVT - SPT;

H0 = H;
% Work out ground state for toric point.
[gs,en] = eigs(H,1,'sa'); gs(abs(gs)<1e-15) = 0;
% Save toric point values separately.
save(['TC N ' num2str(Nv) ' exact ground state.mat'],'basis','H','gs','en');
% Z magnetic field contributions.
MagZ = sparse(1,1,0,2^Nv,2^Nv,2^Nv);
for n = 1:Nv
    MagZ = MagZ + Sz{n}/8;
end
Bz = [0.1 0.05 0.02 0.01];
% X magnetic field contributions.
MagX = sparse(1,1,0,2^Nv,2^Nv,2^Nv);
for n = 1:Nv
    MagX = MagX + Sx{n}/8;
end
Bx = [0.1 0.05 0.02 0.01];

for b = 1:numel(Bz)
    H = H0 + Bz(b)*MagZ;
    [gs,en] = eigs(H,1,'sa');
    save(['TC Bz ' num2str(Bz(b)) ' N ' num2str(Nv) ' exact ground state.mat'],'basis','H','gs','en');
end
for b = 1:numel(Bx)
    H = H0 + Bx(b)*MagX;
    [gs,en] = eigs(H,1,'sa');
    save(['TC Bx ' num2str(Bx(b)) ' N ' num2str(Nv) ' exact ground state.mat'],'basis','H','gs','en');
end