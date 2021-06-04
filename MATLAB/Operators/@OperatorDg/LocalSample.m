% --- Single configuration correlation profile sampling function ---

function [CorrSamp] = LocalSample(Operator,Cfg,EnLoc,dLogp,~)
% This function evaluates the expectation values of an Operator that is
% diagonal in the configuration basis.

CorrSamp = Operator.CfgVal(Operator.Hilbert,Cfg,EnLoc,dLogp,Operator.Graph.Bonds); 

end