% --- Density correlation configuration value function ---

function [CfgVal] = StatSF_Bose_CfgVal(Hilbert,Cfg,~,~,Bonds)
% Given a bosonic configuration Cfg this function computes the value of the
% static structure factor Nq = <Cfg|n{q}n{-q}|Cfg>. Generalisation to more
% than one dimension may require adjustment.
Dim = size(Bonds,2);
[Cfg_vec] = Hilbert.FullCfg(Cfg);
CfgQ_vec = fft(Cfg_vec);
CfgVal = abs(CfgQ_vec).^2;
end