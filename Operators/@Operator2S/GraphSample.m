% --- Two site correlation profile sampling function ---

function [CorrSamp] = GraphSample(Operator,Cfg,~,~,Ansatz)
% This function evaluates the expectation values of a two-site Operator
% with matrix elements detailed in Operator. The Graph version samples over
% all possible Bonds in Graph, not just primary Bonds.

% Make local copies to reduce notation in code below:
BondMap = Operator.Graph.BondMap; CorrSamp = zeros(size(BondMap));

for n = 1:numel(BondMap)
    % Feed the CorrMatEls different elements of BondMap which detail
    % different separations between sites.
    TempGraph.Bonds = BondMap{n};
    [Diff,CorrMatEls] = Operator.CorrMatEls(Operator.Hilbert,Cfg,TempGraph);
    for k = 1:numel(CorrMatEls)
        % Compute ratios with which the matrix elements are 'weighted'.
        PsiRatio = Ansatz.PsiRatio(Ansatz,Diff(k));
        CorrSamp(n) = CorrSamp(n) + PsiRatio * CorrMatEls(k);
    end
end
end