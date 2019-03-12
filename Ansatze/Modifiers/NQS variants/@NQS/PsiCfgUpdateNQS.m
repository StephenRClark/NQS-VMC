% --- General NQS wave function update function ---

function [NQSObj] = PsiCfgUpdateNQS(NQSObj,ThetaP)
% This function updates the intermediate configuration state information
% (effective angles) retained in the NQS ansatz.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP.
% ---------------------------------

% Just overwrite the effective angle information computed earlier.
NQSObj.Theta = ThetaP;
