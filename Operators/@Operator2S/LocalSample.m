% --- Two site correlation profile sampling function ---

function [CorrSamp] = LocalSample(OperatorObj,Cfg,~,~,AnsatzTest)
% This function evaluates the expectation values of a one-site Operator
% with matrix elements detailed in Operator. As single site operators do
% not involve Graph, the Graph and Local Sample functions are the same.

CorrSamp = 0;

[Diff,CorrMatEls] = OperatorObj.CorrMatEls(OperatorObj.Hilbert,Cfg,OperatorObj.Graph);
for k = 1:numel(CorrMatEls)
    % Compute ratios with which the matrix elements are 'weighted'.
    PsiRatio = AnsatzTest.PsiRatio(AnsatzTest,Diff(k));
    CorrSamp = CorrSamp + PsiRatio * CorrMatEls(k);
end

end