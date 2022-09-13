% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSSXTI(NQSObj)
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
% Properties added with translation invariance:
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------
% Format for Params vector is a vertically concatenated stack:
% - (1 x 1) for a.
% - (1 x 1) for A.
% - (Alpha x 1) for b.
% - (Alpha x 1) for B.
% - (Alpha*Nv x 1) for W.
% - (Alpha*Nv x 1) for X.
% ---------------------------------

Nv = NQSObj.Nv; Alpha = NQSObj.Alpha;

Params = zeros(NQSObj.Np,1);

Params(1) = NQSObj.ati;
Params(2) = NQSObj.Ati;
Params((1:Alpha)+2) = NQSObj.bti;
Params((1:Alpha)+2+Alpha) = NQSObj.Bti;
Params((1:(Alpha*Nv))+2+2*Alpha) = reshape(NQSObj.Wv,Alpha*Nv,1);
Params((1:(Alpha*Nv))+2+Alpha*(2+Nv)) = reshape(NQSObj.Xv,Alpha*Nv,1);
end