function [Cfg_vec] = FullBoseCfg(Cfg) 
% Extracts the configuration structure for a spinless boson configuration state
% into a vector of occupation numbers.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-0 bosons here.
% - Cfg.N = total number of sites in the system.
% - Cfg.Nb = total number of bosons in the system.
% - Cfg.occ = (N x 1) vector - boson occupation numbers by site.
% - Cfg.Nmax = maximum number of bosons on a single site.
% ---------------------------------

% Basic field call. Exists as a function for compatibility with other
% general configuration vector representation functions.
Cfg_vec = Cfg.occ;