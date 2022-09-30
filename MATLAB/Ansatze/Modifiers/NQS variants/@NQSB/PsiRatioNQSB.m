% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSB(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQSB Modifier:
% - NQSB.Nv = number of "visible" units.
% - NQSB.Nh = number of "hidden" units.
% - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSB.a = (Nv x 1) vector - visible site bias.
% - NQSB.av = (Nsl x 1) vector - visible bias parameters.
% - NQSB.A = (Nv x 1) vector - visible site square bias.
% - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSB.b = (Nh x 1) vector - hidden site bias.
% - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSB.B = (Nh x 1) vector- hidden site square bias.
% - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSB.HDim = dimension of the hidden units.
% - NQSB.Theta = (Nh x 1) vector - effective angles.
% - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - vector of new effective angles ThetaP.
% Update.NsqVec - vector of new squared visible occupancies.
% ---------------------------------

Ratio = exp(sum(Diff.val.'.*NQSObj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
Theta_shift = zeros(NQSObj.Nh,1); % Initialise effective angle shift.
Nsq_shift = zeros(NQSObj.Nv,1);
% Only loop over the sites where there are differences:
for i=1:Diff.num
    Theta_shift = Theta_shift + Diff.val(i)*NQSObj.W(:,Diff.pos(i));
    Nsq_shift(Diff.pos(i)) = 2*sqrt(NQSObj.NsqVec(Diff.pos(i)))*Diff.val(i) + (Diff.val(i)^2);
end
NsqP = NQSObj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*NQSObj.A(Diff.pos))); % Compute visible square bias contribution.
ThetaP = NQSObj.Theta + Theta_shift; % Update the effective angle for the proposed configuration.
Ratio = Ratio * prod(NHTrace(ThetaP,NQSObj.B,NQSObj.HDim) ./ ...
    NHTrace(NQSObj.Theta,NQSObj.B,NQSObj.HDim)); % Compute full ratio.
% Collect new configuration information into Update.
Update.Theta = ThetaP; Update.NsqVec = NsqP;
end