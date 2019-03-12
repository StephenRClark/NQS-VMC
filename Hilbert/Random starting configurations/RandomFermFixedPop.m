function Cfg = RandomFermFixedPop(Params)
% Creates a random configuration state for a system of N_up + N_dn
% fermions.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

Cfg.N = Params.N;
Cfg.up = sort(randperm(Params.N,Params.N_up));
Cfg.dn = sort(randperm(Params.N,Params.N_dn));

end
