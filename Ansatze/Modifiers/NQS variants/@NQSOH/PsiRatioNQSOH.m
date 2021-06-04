% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSOH(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP and
% new one-hot vector OHVecP.
% ---------------------------------

Nv = NQSObj.Nv; VDim = NQSObj.VDim; OHVec = NQSObj.OHVec;
dV = NQSObj.VList(2) - NQSObj.VList(1); % Needed for half-integer spin cases.
% Convert difference in configuration to difference in one-hot vector.
dOH = zeros(VDim*Nv,1);
for d = 1:Diff.num
    SegInds = (1:VDim) + VDim*(Diff.pos(d)-1); Ind0 = sum((1:VDim).'.*OHVec(SegInds));
    Ind1 = 1+mod(Ind0-1 + (Diff.val(d)/dV),VDim) +  VDim*(Diff.pos(d)-1);
    dOH(SegInds) = -OHVec(SegInds); dOH(Ind1) = 1;
end
ThetaP = NQSObj.Theta + NQSObj.W*dOH; OHVecP = OHVec + dOH;
Ratio = exp(sum(NQSObj.a .* dOH)) * prod(cosh(ThetaP)./cosh(Theta));
Update.OHVec = OHVecP; Update.Theta = ThetaP;
end