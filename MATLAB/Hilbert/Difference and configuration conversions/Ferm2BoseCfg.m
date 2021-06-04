function [Cfg] = Ferm2BoseCfg(Cfg) 
% Converts the configuration structure for a fermionic configuration state
% into a vector of occupation numbers.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------
% Bosonic configurations instead have Cfg_vec stored as Cfg.occ.

Cfg_vec = zeros(Cfg.N,1);
for i=1:numel(Cfg.up)
  Cfg_vec(Cfg.up(i)) = Cfg_vec(Cfg.up(i))+1; 
end
for i=1:numel(Cfg.dn)
  Cfg_vec(Cfg.dn(i)) = Cfg_vec(Cfg.dn(i))+1; 
end

Cfg.occ = Cfg_vec;