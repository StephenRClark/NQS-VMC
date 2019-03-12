% --- Single configuration correlation profile sampling function ---

function [CorrSamp] = GraphSample(Operator,Cfg,EnLoc,dLogp,~)
% This function evaluates the expectation values of an Operator that is
% diagonal in the configuration basis. The Graph version samples over
% all possible Bonds in Graph, not just primary Bonds.

% Make local copies to reduce notation in code below:
BondMap = Operator.Graph.BondMap;

CfgVal = Operator.CfgVal(Operator.Hilbert,Cfg,EnLoc,dLogp,BondMap{1});

CorrSamp = zeros(numel(CfgVal),numel(BondMap)); CorrSamp(:,1) = reshape(CfgVal,numel(CfgVal),1);

for n = 2:numel(BondMap)
    CorrSamp(:,n) = reshape(Operator.CfgVal(Operator.Hilbert,Cfg,EnLoc,dLogp,BondMap{n}),numel(CfgVal),1);
end