% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQS(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
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

Nv = NQSObj.Nv; Nh = NQSObj.Nh;

Params = zeros(NQSObj.Np,1);

Params(1:Nv) = NQSObj.a;
Params((1:Nh)+Nv) = NQSObj.b;
Params((1:(Nh*Nv))+Nh+Nv) = reshape(NQSObj.W,Nh*Nv,1);
end