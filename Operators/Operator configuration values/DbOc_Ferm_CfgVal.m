% --- Doublon number correlation configuration value function ---

function [CfgVal] = DbOc_Ferm_CfgVal(Hilbert,Cfg,~,~,~)
% Given a fermionic configuration Cfg this function computes the matrix
% element <Cfg|D{i}|Cfg>.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.

CfgVal = sum(prod(Cfg_vec,2)); % Compute matrix element.