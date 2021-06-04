% --- Single site correlation profile sampling function ---

function [CorrSamp] = LocalSample(Operator,Cfg,EnLoc,dLogp,AnsatzObj)
% This function evaluates the expectation values of a one-site Operator
% with matrix elements detailed in Operator.

CorrSamp = 0;

% Generate the configuration changes and the Differences linking them to Cfg.
[Diff, CorrMatEls] = Operator.CorrMatEls(Cfg);
for k = 1:numel(CorrMatEls)
    % Compute ratios with which the matrix elements are 'weighted'.
    PsiRatio = AnsatzObj.PsiRatio(Diff(k));
    if isnan(PsiRatio) || isinf(PsiRatio)
        PsiRatio = 0;
    end
    CorrSamp = CorrSamp + PsiRatio * CorrMatEls(k);
end
end