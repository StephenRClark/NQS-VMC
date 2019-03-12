function [Cfg_vec] = FullFermCfg(Cfg)
% Converts the configuration structure for a spin-1/2 fermionic
% configuration state into a vector.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed fermionic here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------
% Cfg_vec is represented as two vectors of up and down occupations.

Cfg_vec = zeros(Cfg.N,2);
for i=1:numel(Cfg.up)
    Cfg_vec(Cfg.up(i),1) = 1;
end
for i=1:numel(Cfg.dn)
    Cfg_vec(Cfg.dn(i),2) = 1;
end