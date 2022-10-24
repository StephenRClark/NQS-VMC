% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSP(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
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
% - (Nsl x VOrder) x 1 for d/da. Group by Sublattice > Visible order
% > [sl,vo], [sl, vo+1] ... [sl+1, vo]
% - (Alpha x HOrder) x 1 for d/db. Group by Alpha > Hidden order
% > [al, ho], [al, ho+1] ... [al+1, ho]
% - (Alpha x Nv) x (HOrder x VOrder) for d/dW. Group by Alpha > Position > Hidden order > Visible order
% > [al,v,ho,vo], [al,v,ho,vo+1] ... [al,v,ho+1,vo] ... [al,v+1,ho,vo] ... [al+1,v,ho,vo]
% ---------------------------------

Params = zeros(NQSObj.Np,1);
W_shift = permute(NQSObj.W,[4 3 2 1]);
a_vec = NQSObj.a.'; b_vec = NQSObj.b.'; 
Params = Params + [a_vec(:); b_vec(:), W_shift(:)]; % Will throw error if not same length.
end