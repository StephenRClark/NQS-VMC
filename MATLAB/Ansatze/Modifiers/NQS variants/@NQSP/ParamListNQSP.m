% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSP(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
% interactions. Max value VDim-1.
% - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
% interactions. Max value HDim-1.
% - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
% HOrder) + (Nv x VOrder)(Nh x HOrder)
% - NQS.a = (Nv x VOrder) matrix - visible site biases.
% - NQS.b = (Nh x HOrder) matrix - hidden site bias.
% - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.VisVec = (Nv x 1) vector - visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x VOrder) x 1 for d/da.
% - (Nh x HOrder) x 1 for d/db.
% - (Nh x Nv) x (HOrder x VOrder) for d/dW.
% ---------------------------------

Nv = NQSObj.Nv; VOrder = NQSObj.VOrder;
Nh = NQSObj.Nh; HOrder = NQSObj.HOrder;

Params = zeros(NQSObj.Np,1);

Params(1:Nv*VOrder) = NQSObj.a(:);
Params((1:Nh*HOrder)+Nv*VOrder) = NQSObj.b(:);
Params((1:(Nh*Nv))+Nh*HOrder+Nv*VOrder) = NQSObj.W(:);
end