% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSS1(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

Nv = NQSObj.Nv; Nh = NQSObj.Nh;

Params = zeros(NQSObj.Np,1);

Params(1:Nv) = NQSObj.a;
Params((1:Nv)+Nv) = NQSObj.A;
Params((1:Nh)+2*Nv) = NQSObj.b;
Params((1:(Nh*Nv))+2*Nv+Nh) = reshape(NQSObj.w,Nh*Nv,1);
Params((1:(Nh*Nv))+2*Nv+Nh+Nv*Nh) = reshape(NQSObj.W,Nh*Nv,1);
end