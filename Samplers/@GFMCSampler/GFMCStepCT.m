% --- Continuous time Green's Function Monte Carlo step function ---

function [EnLoc,OpLoc,AnsatzObj,Cfg,Weight] = GFMCStepCT(HamiltonianObj,AnsatzObj,Cfg,T,Weight,Operators)

% Assumes AnsatzObj has been appropriately prepared with input Cfg.
OpLoc = cell(1,numel(Operators)); TRem = T; AnsatzObj = AnsatzObj.PrepPsi(Cfg);
while TRem > 0
    [Diff,HMatEls] = HamiltonianObj.HamMatEls(Cfg);
    [CfgP] = AnsatzObj.Hilbert.Diff2Cfg(Diff,Cfg);
    Ratios = zeros(size(HMatEls,1),1); Updates = cell(size(HMatEls,1),1);
    DiagInds = zeros(size(HMatEls,1),1); EnLoc = 0;
    for d = 1:numel(HMatEls)
        [Ratios(d), Updates{d}] = AnsatzObj.PsiRatio(Diff(d));
        HMatEls(d) = HMatEls(d) * Ratios(d);
        EnLoc = EnLoc + HMatEls(d);
        if Diff(d).num == 0 % Record which Diffs correspond to diagonal moves.
            DiagInds(d) = 1;
        end
    end
    GDiag = sum(HMatEls(DiagInds==1));
    R = rand; KT = log(1-R)/(EnLoc - GDiag); DT = min(TRem,KT);
    TRem = TRem - DT; Weight = Weight * exp(-DT*EnLoc);
    if TRem > 0 % Remaining time suggests a transition occurs.
        Px = -HMatEls(DiagInds==0); UInds = 1:size(HMatEls,1); UInds = UInds(DiagInds==0);
        PxVec = cumsum(Px) / sum(Px); R = rand; Ind = sum((PxVec-R)<0)+1;
        Cfg = CfgP(UInds(Ind)); AnsatzObj = PsiCfgUpdate(AnsatzObj,Updates{UInds(Ind)});
    end
end
% Calculate local values of energy and other operators.
for o = 1:numel(Operators)
    OpLoc{o} = Operators{o}.GraphSample(Cfg,EnLoc,0,AnsatzObj);
end
end