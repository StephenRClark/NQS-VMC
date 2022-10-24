% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSS1(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQSS1.Nv = number of "visible" spins.
% - NQSS1.Nh = number of "hidden" spins.
% - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
% - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSS1.a = (Nv x 1) vector - visible site bias.
% - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
% - NQSS1.A = (Nv x 1) vector - visible site square bias.
% - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSS1.b = (Nh x 1) vector - hidden site bias.
% - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
% - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
% - NQSS1.Theta = (Nh x 1) vector - effective angles.
% - NQSS1.VisVec = (Nv x 1) vector - visible occupancies vector.
% - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - vector of new effective angles ThetaP.
% Update.VisVec - vector of new visible occupancies.
% Update.NsqVec - vector of new squared visible occupancies.
% ---------------------------------

Ratio = exp(sum(Diff.val.'.*NQSObj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
Theta_shift = zeros(NQSObj.Nh,1); % Initialise effective angle shift.
Nsq_shift = zeros(NQSObj.Nv,1);
% Only loop over the sites where there are differences:
for i=1:Diff.num
    Theta_shift = Theta_shift + Diff.val(i)*NQSObj.w(:,Diff.pos(i));
    Nsq_shift(Diff.pos(i)) = -2*NQSObj.NsqVec(Diff.pos(i)) + (Diff.val(i)^2);
    Theta_shift = Theta_shift + Nsq_shift(Diff.pos(i))*NQSObj.W(:,Diff.pos(i));
end
NsqP = NQSObj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*NQSObj.A(Diff.pos))); % Compute visible square bias contribution.
ThetaP = NQSObj.Theta + Theta_shift; % Update the effective angle for the proposed configuration.
Ratio = Ratio * prod(cosh(ThetaP) ./ cosh(NQSObj.Theta)); % Compute full ratio.
% Collect new configuration information into Update.
Update.Theta = ThetaP; Update.NsqVec = NsqP;
end