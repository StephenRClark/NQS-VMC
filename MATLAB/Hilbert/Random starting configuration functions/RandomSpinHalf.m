function Cfg = RandomSpinHalf(CParams)
% Creates a random configuration state for a system of N spin-1/2's.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

N = CParams.N; Cfg.N = N;

Cfg_vec = zeros(N,1);
for n = 1:N
    Cfg_vec(n) = (-1)^(mod(randi(2),2));
end
Cfg.up = find(Cfg_vec>0).'; Cfg.dn = find(Cfg_vec<0).';
end