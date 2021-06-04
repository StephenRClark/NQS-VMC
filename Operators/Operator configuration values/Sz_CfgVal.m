% --- Single site spin correlation configuration value function ---

function [CfgVal] = Sz_CfgVal(Hilbert,Cfg,~,~,~)
% Given a spin-1/2 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|Sz{i}|CfgP>.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.

% This function gives total Sz value - for a profile, use DenProf.
CfgVal = sum(Cfg_vec);