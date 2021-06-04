function Cfg = RandomUnarySpin(CParams)
% Creates a random configuration state for a system of N spin-1/2's.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% - Cfg.N0 = original number of sites.
% - Cfg.d = original site dimension.
% ---------------------------------

N = CParams.N; N0 = CParams.N0; d = CParams.d; 
Cfg.N = N; Cfg.N0 = N0; Cfg.d = d;

Cfg_mat = -ones(d,N0);
for n = 1:N0
    v0 = randi(d);
    Cfg_mat(v0,n) = 1;
end
Cfg_vec = reshape(Cfg_mat,N,1);
Cfg.up = find(Cfg_vec>0).'; Cfg.dn = find(Cfg_vec<0).';
end