% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSSH(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nh + 2*Alpha + 2.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.B = (Nh x 1) vector - hidden site square bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
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
Ratio = Ratio * prod(SHTrace(ThetaP,NQSObj.B,NQSObj.HDim) ./ ...
    SHTrace(NQSObj.Theta,NQSObj.B,NQSObj.HDim)); % Compute full ratio.
% Collect new configuration information into Update.
Update.Theta = ThetaP; Update.NsqVec = NsqP;