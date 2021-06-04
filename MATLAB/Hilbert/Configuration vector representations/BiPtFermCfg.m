function [Cfg_vec] = BiPtFermCfg(Cfg)
% Converts the configuration structure for a spin-1/2 fermionic
% configuration state into a vector.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed fermionic here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% - Cfg.db = (Ndb x 1) vector of sites with double occupations.
% - Cfg.mt - (Nmt x 1) vector of sites with no particles.
% ---------------------------------
% Where FullFermCfg creates a vector listing empty, doubly occupied and
% single fermion sites, this function creates a vector of 
% [ ---n_up--- ---n_dn--- ] occupations, with + meaning occupied. Used for 
% bipartite NQS structures.

Cfg_vec = zeros(2*Cfg.N,1);
for i=1:numel(Cfg.up)
    Cfg_vec(Cfg.up(i)) = 1;
end
for i=1:numel(Cfg.dn)
    Cfg_vec(Cfg.dn(i)+Cfg.N) = 1;
end