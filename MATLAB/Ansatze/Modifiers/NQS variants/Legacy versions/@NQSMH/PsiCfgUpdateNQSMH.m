% --- General NQS wave function update function ---

function [NQSObj] = PsiCfgUpdateNQSMH(NQSObj,Update)
% This function updates the intermediate configuration state information
% (effective angles) retained in the NQS ansatz.
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

% Just overwrite the configuration information computed earlier.
NQSObj.ThetaH = Update.ThetaH; NQSObj.ThetaM = Update.ThetaM;
NQSObj.NsqVec = Update.NsqVec; NQSObj.Hv = Update.Hv; NQSObj.Mv = Update.Mv;
