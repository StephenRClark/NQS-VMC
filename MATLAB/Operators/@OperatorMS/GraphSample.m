% --- Multi site correlation profile sampling function ---

function [CorrSamp] = GraphSample(OperatorObj,Cfg,~,~,AnsatzObj)
% This function evaluates the expectation values of a multi-site Operator
% with matrix elements detailed in Operator. 

% The multi-site GraphSample will have separate entries for each cluster
% specified in the Operator Graph's ExtraLabels list.

CorrSamp = cell(size(OperatorObj.Graph.ExtraLists{OperatorObj.ListIndex},1),1);

% Generate the configuration changes and the Differences linking them to Cfg.
[Diff, CorrMatEls] = OperatorObj.CorrMatEls(Cfg);
for k = 1:numel(CorrMatEls)
    % Compute ratios with which the matrix elements are 'weighted'.
    PsiRatio = AnsatzObj.PsiRatio(Diff(k));
    if isnan(PsiRatio) || isinf(PsiRatio)
        PsiRatio = 0;
    end
    % Each Diff will be labelled by an index indicating which site grouping
    % it is relevant for.
    CorrSamp{Diff(k).index} = PsiRatio * CorrMatEls(k);
end
end