% --- Single configuration correlation profile sampling function ---

function [Diff, OpMatEls] = CorrMatEls(OperatorObj,HilbertObj,Cfg,GraphObj)
% This function evaluates the expectation values of an Operator that is
% diagonal in the configuration basis, and also outputs a zero Difference
% struct - used when CorrMatEls is requested of a diagonal operator.

% Generate a zero Difference struct.
Diff.num = 0;
Diff.pos = 1;
Diff.val = 0;
Diff.type = 0;
Diff.sign = 1;

OpMatEls = OperatorObj.CfgVal(HilbertObj,Cfg,0,0,GraphObj.Bonds); 

end