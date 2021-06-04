% --- Two site correlation profile sampling function ---

function [CorrSamp] = GraphSample(OperatorObj,Cfg,~,~,AnsatzObj)
% This function evaluates the expectation values of a two-site Operator
% with matrix elements detailed in Operator. The Graph version samples over
% all possible Bonds in Graph, not just primary Bonds.

% Make local copies to reduce notation in code below:
Graph = OperatorObj.Graph; BondMap = Graph.BondMap; CorrSamp = cell(Cfg.N);

for n = 1:numel(BondMap)
    % Feed the CorrMatEls different elements of BondMap which detail
    % different separations between sites.
    TempGraph = Graph; TempGraph.Bonds = BondMap{n}; OperatorObj.Graph = TempGraph;
    [Diff,CorrMatEls] = OperatorObj.CorrMatEls(Cfg);
    if isempty(CorrMatEls) == 0
        CorrMatElsRef = CorrMatEls(1,:);
    end
    for k = 1:size(CorrMatEls,1)
        % Compute ratios with which the matrix elements are 'weighted'.
        PsiRatio = AnsatzObj.PsiRatio(Diff(k));
        if isnan(PsiRatio) || isinf(PsiRatio)
            PsiRatio = 0;
        end
        if Diff(k).num == 0 % Sometimes occurs if matrix element function allows for same site operation.
            P1 = Diff(k).pos; P2 = Diff(k).pos;
        else
            P1 = Diff(k).pos(1); P2 = Diff(k).pos(2);
        end
        if isempty(CorrSamp{P1,P2})
            CorrSamp{P1,P2} = PsiRatio * CorrMatEls(k,:);
        else
            CorrSamp{P1,P2} = CorrSamp{P1,P2} + PsiRatio * CorrMatEls(k,:);
        end
    end
end
for n = 1:Cfg.N
    for m = 1:Cfg.N
        if isempty(CorrSamp{n,m})
            CorrSamp{n,m} = 0 * CorrMatElsRef;
        end
    end
end
CorrSamp = cell2mat(CorrSamp);
end