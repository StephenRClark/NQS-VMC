function Cfg = RandomSpin1NoSector(CParams)
% Creates a random configuration state for a system of N spin-1s 
% with no restriction on total spin projection.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.SzT = total configuration Sz.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

N = CParams.N; Cfg.N = N;
Cfg_vec = randi(3,[N,1]) - 2;
Cfg.up = find(Cfg_vec == 1).'; Cfg.dn = find(Cfg_vec == -1).';
end

