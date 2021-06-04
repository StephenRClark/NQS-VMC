% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSTI(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQS Modifier object with translation invariance:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nh + Alpha + 1.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% Properties added with translation invariance:
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
Nv = NQSObj.Nv; Alpha = NQSObj.Alpha;

Params = zeros(NQSObj.Np,1);

Params(1) = NQSObj.ati;
Params((1:Alpha)+1) = NQSObj.bti;
Params((1:(Alpha*Nv))+1+Alpha) = reshape(NQSObj.Wv.',Alpha*Nv,1);
end