% --- Two site correlation profile sampling function ---

function [CorrSamp] = LocalSample(OperatorObj,Cfg,~,~,AnsatzObj)
% This function evaluates the expectation values of a two-site Operator
% with matrix elements detailed in Operator.

[Diff,CorrMatEls] = OperatorObj.CorrMatEls(Cfg);
if numel(CorrMatEls) ~= 0
    for k = 1:numel(CorrMatEls)
        % Compute ratios with which the matrix elements are 'weighted'.
        PsiRatio = AnsatzObj.PsiRatio(Diff(k));
        if isnan(PsiRatio) || isinf(PsiRatio)
            PsiRatio = 0;
        end
        if k == 1
            CorrSamp = PsiRatio * CorrMatEls(k,:);
        else
            CorrSamp = CorrSamp + PsiRatio * CorrMatEls(k,:);
        end
    end
else
    CorrSamp = 0;
end

end