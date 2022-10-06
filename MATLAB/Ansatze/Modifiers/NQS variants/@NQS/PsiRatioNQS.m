% --- General NQS wave function amplitude ratio function ---

function [Ratio,ThetaP] = PsiRatioNQS(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.Alpha = number of unique coupling sets or "hidden unit density".
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.av = (Nsl x 1) vector - visible bias parameters.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP.
% ---------------------------------

Ratio = exp(sum(Diff.val.'.*NQSObj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
Theta_shift = zeros(NQSObj.Nh,1); % Initialise effective angle shift.
% Only loop over the sites where there are differences:
for i=1:Diff.num
  Theta_shift = Theta_shift + Diff.val(i)*NQSObj.W(:,Diff.pos(i));
end
ThetaP = NQSObj.Theta + Theta_shift; % Update the effective angle for the proposed configuration.
Ratio = Ratio * prod(cosh(ThetaP)./cosh(NQSObj.Theta)); % Compute full ratio.
% Return ThetaP information in case proposed move is accepted for speedy ansatz configuration update.
end