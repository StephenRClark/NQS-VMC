% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSS1TI(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Alpha + 2 + Alpha.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% Properties added with translation invariance:
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------
Nv = NQSObj.Nv; Alpha = NQSObj.Alpha;

Params = zeros(NQSObj.Np,1);

Params(1) = NQSObj.ati;
Params(2) = NQSObj.Ati;
Params((1:Alpha)+2) = NQSObj.bti;
Params((1:(Alpha*Nv))+2+Alpha) = reshape(NQSObj.wv,Alpha*Nv,1);
Params((1:(Alpha*Nv))+2+Alpha*(1+Nv)) = reshape(NQSObj.Wv,Alpha*Nv,1);
end