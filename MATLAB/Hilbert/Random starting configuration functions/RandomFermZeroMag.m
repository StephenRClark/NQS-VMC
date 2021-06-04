function Cfg = RandomFermZeroMag(Params)
% Creates a random configuration state for a system of Nf fermions
% with zero magnetisation.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

if mod(Params.Nf,2)==0
    Cfg.N = Params.N;
    Cfg.up = sort(randperm(Params.N,Params.Nf/2));
    Cfg.dn = sort(randperm(Params.N,Params.Nf/2));
else
    error('RandomFermZeroMag(N,Nf): Even number of fermions Nf required.')
end
