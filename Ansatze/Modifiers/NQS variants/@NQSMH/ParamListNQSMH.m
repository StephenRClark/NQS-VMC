% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSMH(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
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

Nv = NQSObj.Nv; Nh = NQSObj.Nh;

Params = zeros(NQSObj.Np,1);

Params(1:Nv) = NQSObj.a;
Params((1:Nv)+Nv) = NQSObj.A;
Params((1:Nh)+2*Nv) = NQSObj.BH;
Params((1:Nh)+2*Nv+Nh) = NQSObj.BM;
Params((1:(Nh*Nv))+2*(Nv+Nh)) = reshape(NQSObj.W,Nh*Nv,1);
Params((1:(Nh*Nv))+2*(Nv+Nh)+Nh*Nv) = reshape(NQSObj.X,Nh*Nv,1);
end