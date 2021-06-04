% --- Constant onfiguration value function ---

function [CfgVal] = Const_CfgVal(~,Cfg,~,~,~)
% A constant operator for adding offsets to energy terms. Outputs the
% number of sites in the configuration.

CfgVal = Cfg.N;