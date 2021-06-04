% --- Two site density correlation configuration value function ---

function [CfgVal] = VarN_Bose_CfgVal(Hilbert,Cfg,~,~,~)
% Given a bosonic configuration Cfg this function computes the matrix
% element <Cfg|N^2|Cfg> - <Cfg|N|Cfg>^2.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Nb = Cfg.Nb; % Number of bosons in the configuration.

% This operator is diagonal in the boson number basis.
CfgVal = sum((Cfg_vec - (Nb/N)).^2)/N;