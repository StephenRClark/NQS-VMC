% --- Density profile configuration value function ---

function [CfgVal] = DenSum_CfgVal(Hilbert,Cfg,~,~,~)
% Given any configuration Cfg this function computes the configuration
% density profile in terms of the single-site operator characterising the
% configuration. 

% Spin Hilbert - Sz profile.
% Ferm Hilbert - Up and down density profiles.
% Bose Hilbert - Density profile.

Cfg_vec = Hilbert.FullCfg(Cfg);

CfgVal = sum(Cfg_vec(:));