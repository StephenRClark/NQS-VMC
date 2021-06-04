% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSMH(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nv + 2*Nh + 2*Nv.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
% - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------
% Format for Update is a struct with five fields:
% Update.ThetaH - vector of new effective angles ThetaHP.
% Update.ThetaM - vector of new effective angles ThetaMP.
% Update.NsqVec - vector of new squared visible occupancies.
% Update.Hv - vector of new holon operator values HvP.
% Update.Mv - vector of new multiplon operator values MvP.
% ---------------------------------

Ratio = exp(sum(Diff.val.'.*NQSObj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
ThetaH_shift = zeros(NQSObj.Nh,1); ThetaM_shift = zeros(NQSObj.Nh,1); % Initialise effective angle shifts.
Nsq_shift = zeros(NQSObj.Nv,1); Hv = NQSObj.Hv; Mv = NQSObj.Mv; HDim = NQSObj.HDim;
dHv = zeros(NQSObj.Nv,1); dMv = zeros(NQSObj.Nv,1);
% Only loop over the sites where there are differences:
for i=1:Diff.num
    dHv(Diff.pos(i)) = ((Mv(Diff.pos(i)) + Diff.val(i)) < 0) - (Diff.val(i) > 0)*Hv(Diff.pos(i));
    dMv(Diff.pos(i)) = Diff.val(i) + dHv(Diff.pos(i)); % Holon / multiplon operator differences.
    Nsq_shift(Diff.pos(i)) = 2*sqrt(NQSObj.NsqVec(Diff.pos(i)))*Diff.val(i) + (Diff.val(i)^2);
    ThetaH_shift = ThetaH_shift + dHv(Diff.pos(i))*NQSObj.W(:,Diff.pos(i)) + ...
        dMv(Diff.pos(i))*NQSObj.X(:,Diff.pos(i)); % Change in hidden holon effective angles.
    ThetaM_shift = ThetaM_shift + dMv(Diff.pos(i))*NQSObj.W(:,Diff.pos(i)) + ...
        dHv(Diff.pos(i))*NQSObj.X(:,Diff.pos(i)); % Change in hidden multiplon effective angles.
end
NsqP = NQSObj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
HvP = Hv + dHv; MvP = Mv + dMv; % Update the holon / multiplon operator vectors for the proposed configuration.
Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*NQSObj.A(Diff.pos))); % Compute visible square bias contribution.
ThetaHP = NQSObj.ThetaH + ThetaH_shift; ThetaMP = NQSObj.ThetaM + ThetaM_shift;
% Update the effective angles for the proposed configuration.
Ratio = Ratio * prod(MHTrace(ThetaHP,ThetaMP,HDim) ./ MHTrace(NQSObj.ThetaH,NQSObj.ThetaM,HDim)); % Compute full ratio.
% Collect new configuration information into Update.
Update.ThetaH = ThetaHP; Update.ThetaM = ThetaMP;
Update.NsqVec = NsqP; Update.Hv = HvP; Update.Mv = MvP;