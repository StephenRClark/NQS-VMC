% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSOH(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% ---------------------------------

Nv = NQSObj.Nv; Nh = NQSObj.Nh; VDim = NQSObj.VDim;

Params = zeros(NQSObj.Np,1);

Params(1:(VDim*Nv)) = NQSObj.a;
Params((1:Nh)+VDim*Nv) = NQSObj.b;
Params((1:(Nh*VDim*Nv))+Nh+VDim*Nv) = reshape(NQSObj.W,Nh*VDim*Nv,1);
end