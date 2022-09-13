% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSU(NQSObjNew,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" spins.
% - NQSU.Nh = number of "hidden" spins.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP and
% new unary vector UVecP.
% ---------------------------------

Nv = NQSObjNew.Nv; Nmax = NQSObjNew.VDim-1; UVec = NQSObjNew.UVec;
dV = NQSObjNew.VList(2) - NQSObjNew.VList(1); % Needed for half-integer spin cases.
% Convert difference in configuration to difference in unary vector.
dU = zeros(Nmax*Nv,1);
for d = 1:Diff.num
    SegInds = (1:Nmax) + Nmax*(Diff.pos(d)-1); Ind0 = 1+sum((1:Nmax).'.*UVec(SegInds));
    IndP = 1+mod(Ind0-1 + (Diff.val(d)/dV),Nmax);
    dU(SegInds) = -UVec(SegInds);
    if (IndP-1) > 0
        dU(IndP-1 + Nmax*(Diff.pos(d)-1)) = 1;
    end
end
ThetaP = NQSObjNew.Theta + NQSObjNew.W*dU; UVecP = UVec + dU;
Ratio = exp(sum(NQSObjNew.a .* dU)) * prod(cosh(ThetaP)./cosh(Theta));
Update.UVec = UVecP; Update.Theta = ThetaP;
end