% --- Charge density rigidity correlation configuration value function ---

function [CfgVal] = Rgdy_Ferm_CfgVal(Hilbert,Cfg,~,~,~)
% Given a fermionic configuration Cfg this function computes the matrix
% element <Cfg|exp(2i*pi*jn(j)/N)|Cfg>. This will require modification in
% order to apply to systems in higher dimensions.

N = Cfg.N;

[Cfg_vec] = sum(Hilbert.FullCfg(Cfg),2); % Build the configuration vector for Cfg for convenience below.

CfgVal = exp(sum(2i*pi*(Cfg_vec.*((1:N).'))/N));