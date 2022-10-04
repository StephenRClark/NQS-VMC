% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSM(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.VDim = dimensions of the visible units.
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP, new
% visible occupancies VisVec and new species matrix OMatP.
% ---------------------------------

Nv = NQSObj.Nv; OMat = NQSObj.OMat;
VisVecP = NQSObj.VisVec; OMatP = OMat; ThetaP = NQSObj.Theta;
for d = 1:Diff.num
    VisVecP(Diff.pos(d)) = VisVecP(Diff.pos(d)) + Diff.val(d);
    switch VisVecP(Diff.pos(d))
        case 0
            OMatP(Diff.pos(d),:) = [1 0 0];
        case 1
            OMatP(Diff.pos(d),:) = [0 0 0];
        case 2
            OMatP(Diff.pos(d),:) = [0 1 0];
        otherwise
            OMatP(Diff.pos(d),:) = [0 0 1];
    end
end
dHDM = reshape((OMatP-NQSObj.OMat),3*Nv,1); dHD = dHDM(1:(2*Nv));
ThetaP = ThetaP + [NQSObj.W, NQSObj.X] * dHD;
Ratio = exp(sum(NQSObj.a .* dHDM)) * prod(cosh(ThetaP)./cosh(NQSObj.Theta));
Update.OMat = OMatP; Update.VisVec = VisVecP; Update.Theta = ThetaP;
end