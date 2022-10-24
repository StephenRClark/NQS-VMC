% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSSX(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
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
% Format for Params vector is a vertically concatenated stack:
% - (Nv x 1) for a.
% - (Nv x 1) for A.
% - (Nh x 1) for b.
% - (Nh x 1) for B.
% - (Nh*Nv x 1) for W.
% - (Nh*Nv x 1) for X.
% ---------------------------------
Nv = NQSObj.Nv; Nh = NQSObj.Nh;

Params = zeros(NQSObj.Np,1);

Params(1:Nv) = NQSObj.a;
Params((1:Nv)+Nv) = NQSObj.A;
Params((1:Nh)+2*Nv) = NQSObj.b;
Params((1:Nh)+2*Nv+Nh) = NQSObj.B;
Params((1:(Nh*Nv))+2*Nv+2*Nh) = reshape(NQSObj.W,Nh*Nv,1);
Params((1:(Nh*Nv))+2*Nv+2*Nh+Nv*Nh) = reshape(NQSObj.X,Nh*Nv,1);
end