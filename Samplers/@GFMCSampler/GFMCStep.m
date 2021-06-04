% --- Continuous time Green's Function Monte Carlo step function ---

function [EnLoc,OpLoc,AnsatzObj,Cfg,Weight] = GFMCStep(HamiltonianObj,AnsatzObj,Cfg,Lambda,Weight,Operators)

% Assumes AnsatzObj has been appropriately prepared with input Cfg.
OpLoc = cell(1,numel(Operators));
[Diff,HMatEls] = HamiltonianObj.HamMatEls(Cfg);
[CfgP] = AnsatzObj.Hilbert.Diff2Cfg(Diff,Cfg);
Ratios = zeros(size(HMatEls,1),1); Updates = cell(numel(HMatEls),1);
EnLoc = 0; GDiag = Lambda; DiagInds = zeros(size(HMatEls,1),1);
for d = 1:numel(HMatEls)
    [Ratios(d), Updates{d}] = AnsatzObj.PsiRatio(Diff(d));
    EnLoc = EnLoc + (HMatEls(d) * Ratios(d));
    HMatEls(d) = - HMatEls(d) * Ratios(d);    
    if Diff(d).num == 0 % Record which Diffs correspond to diagonal moves.
        GDiag = GDiag + HMatEls(d); DiagInds(d) = 1;
        HMatEls(d) = HMatEls(d) + Lambda;
    end    
end

if sum(DiagInds) == 0 % No diagonal matrix element of H.
    HMatEls = [HMatEls; Lambda]; DiagInds = [DiagInds; 1];
end
Bx = sum(HMatEls); Px = cumsum(HMatEls/Bx); 
Weight = Weight * Bx;

% Calculate local values of energy and other operators.
% EnLoc = HamiltonianObj.EnergySample(Cfg,AnsatzObj);
for o = 1:numel(Operators)
    OpLoc{o} = Operators{o}.GraphSample(Cfg,EnLoc,0,AnsatzObj);
end

R = rand; Ind = sum(Px<R)+1; 
if DiagInds(Ind) == 0 % Non-diagonal move chosen.
    Cfg = CfgP(Ind); AnsatzObj = PsiCfgUpdate(AnsatzObj,Updates{Ind});
end

end