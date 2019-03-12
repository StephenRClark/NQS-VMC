% --- Single site correlation profile sampling function ---

function [CorrSamp] = LocalSample(Operator,Cfg,EnLoc,dLogp,Ansatz)
% This function evaluates the expectation values of a one-site Operator
% with matrix elements detailed in Operator.

CorrSamp = 0;

% Generate the configuration changes and the Differences linking them to Cfg.
[Diff, CorrMatEls] = Operator.CorrMatEls(Operator.Hilbert,Cfg);
for k = 1:numel(CorrMatEls)
    % Compute ratios with which the matrix elements are 'weighted'.
    PsiRatio = Ansatz.PsiRatio(Ansatz,Diff(k));
    CorrSamp = CorrSamp + PsiRatio * CorrMatEls(k);
end
end