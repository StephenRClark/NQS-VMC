% --- Single site correlation profile sampling function ---

function [CorrSamp] = GraphSample(Operator,Cfg,~,~,Ansatz)
% This function evaluates the expectation values of a one-site Operator
% with matrix elements detailed in Operator. As single site operators do
% not involve Graph, the Graph and Local Sample functions are the same.

CorrSamp = 0;

% Generate the configuration changes and the Differences linking them to Cfg.
[Diff, CorrMatEls] = Operator.CorrMatEls(Operator.Hilbert,Cfg);
for k = 1:numel(CorrMatEls)
    % Compute ratios with which the matrix elements are 'weighted'.
    PsiRatio = Ansatz.PsiRatio(Ansatz,Diff(k));
    CorrSamp = CorrSamp + PsiRatio * CorrMatEls(k);
end
end