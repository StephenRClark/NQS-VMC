Nt = 4; Nv = Nt*3;

Triplet = zeros(Nt,3);
LeftQuad = zeros(Nt,4);
RightQuad = zeros(Nt,4);
for t = 1:Nt
    Triplet(t,:) = [1 2 3] + (t-1)*3;
    LeftQuad(t,:) = mod([-1 0 1 2] + (t-1)*3,Nv) + 1;
    RightQuad(t,:) = mod([0 1 2 3] + (t-1)*3,Nv) + 1;
end

basis = zeros(2^Nv,Nv); W = 2^Nv; uinds = zeros(W,1);
for k=1:2^Nv
    v = 3-2*ind2dig(k,2*ones(1,Nv));
    basis(k,:) = v; % Express basis as sz = {+1,-1}
    if sum(v) == -Nt
        % Check if basis follows unary encoding using v
        vm = reshape(v,3,Nt); vs = sum(vm,1);
        if sum(vs==-1)==Nt
            uinds(k) = 1;
        end
    end
end
uinds = uinds>0; basis_u = basis(uinds,:);

[Sx,Sy,Sz,~,~] = spin_half_operators(Nv);
ID = speye(2^Nv);

XYX = sparse(1,1,0,W,W,W);
LXZZX = sparse(1,1,0,W,W,W);
RXZZX = sparse(1,1,0,W,W,W);
for t = 1:Nt
    xyx = Sx{Triplet(t,1)}*Sy{Triplet(t,2)}*Sx{Triplet(t,3)}*8; 
    lxzzx = Sx{LeftQuad(t,1)}*Sz{LeftQuad(t,2)}*Sz{LeftQuad(t,3)}*Sx{LeftQuad(t,4)}*16; 
    rxzzx = Sx{RightQuad(t,1)}*Sz{RightQuad(t,2)}*Sz{RightQuad(t,3)}*Sx{RightQuad(t,4)}*16;
    XYX = XYX + xyx; LXZZX = LXZZX + lxzzx; RXZZX = RXZZX + rxzzx;
end

H = - XYX - LXZZX - RXZZX; 

[gs,en] = eigs(H,1,'sr'); gs = gs/sign(gs(1)); gs(abs(gs)<1e-15) = 0;

gs_u = gs(uinds); gs_u = gs_u/sign(gs_u(1)); gs_u = gs_u/sqrt(sum(abs(gs_u).^2));

save(['GS 1D XYX-LXZZX-RXZZX N ' num2str(Nv) ' exact ground state.mat'],...
    'H','en','basis','basis_u','gs','gs_u','uinds');