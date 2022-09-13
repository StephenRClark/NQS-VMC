% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSOHTI(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim + Alpha + (VDim*Nv * Alpha).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% Properties added with translation invariance:
% - NQS.ati = (VDim x 1) vector - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x VDim*Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------

Nv = NQSObj.Nv; Alpha = NQSObj.Alpha; VDim = NQSObj.VDim;

Params = zeros(NQSObj.Np,1);

Params(1:(VDim)) = NQSObj.ati;
Params((1:Alpha)+VDim) = NQSObj.bti;
Params((1:(Alpha*VDim*Nv))+Alpha+VDim) = reshape(NQSObj.Wv.',Alpha*VDim*Nv,1);
end