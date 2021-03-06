function [Cfg_vec] = FullFermDen(Cfg)
% Converts the configuration structure for a spin-1/2 fermionic
% configuration state into a vector.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed fermionic here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------
% Cfg_vec is represented as a vector of particle occupations, ignoring spin.

Cfg_vec = zeros(Cfg.N,1);
for i=1:numel(Cfg.up)
    Cfg_vec(Cfg.up(i)) = Cfg_vec(Cfg.up(i))+1;
end
for i=1:numel(Cfg.dn)
    Cfg_vec(Cfg.dn(i)) = (Cfg.dn(i))+1;
end