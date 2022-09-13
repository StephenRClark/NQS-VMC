% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSSX(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - vector of new effective angles ThetaP.
% Update.VisVec - vector of new visible occupancies.
% Update.ThetaSq - vector of new effective angles ThetaSqP.
% Update.NsqVec - vector of new squared visible occupancies.
% ---------------------------------

Ratio = exp(sum(Diff.val.'.*NQSObj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
Theta_shift = zeros(NQSObj.Nh,1); ThetaSq_shift = zeros(NQSObj.Nh,1); % Initialise effective angle shift.
Vis_shift = zeros(NQSObj.Nv,1); Nsq_shift = zeros(NQSObj.Nv,1);
% Only loop over the sites where there are differences:
for i=1:Diff.num
    Vis_shift(Diff.pos(i)) = Diff.val(i);
    Theta_shift = Theta_shift + Diff.val(i)*NQSObj.W(:,Diff.pos(i));
    Nsq_shift(Diff.pos(i)) = 2*Diff.val(i)*NQSObj.VisVec(Diff.pos(i)) + Diff.val(i)^2;
    ThetaSq_shift = ThetaSq_shift + Nsq_shift(Diff.pos(i))*NQSObj.X(:,Diff.pos(i));
end
VisP = NQSObj.VisVec + Vis_shift; NsqP = NQSObj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*NQSObj.A(Diff.pos))); % Compute visible square bias contribution.
ThetaP = NQSObj.Theta + Theta_shift; ThetaSqP = NQSObj.ThetaSq + ThetaSq_shift;% Update the effective angle for the proposed configuration.
Ratio = Ratio * prod(SqTrace(ThetaP,ThetaSqP,NQSObj.HVals) ./ ...
    SqTrace(NQSObj.Theta,NQSObj.ThetaSq,NQSObj.HVals)); % Compute full ratio.
% Collect new configuration information into Update.
Update.Theta = ThetaP; Update.ThetaSq = ThetaSqP; Update.VisVec = VisP; Update.NsqVec = NsqP;end