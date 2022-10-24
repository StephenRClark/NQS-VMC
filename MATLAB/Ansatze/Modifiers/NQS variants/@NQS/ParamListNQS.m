% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQS(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.Alpha = number of unique coupling sets or "hidden unit density".
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.av = (Nsl x 1) vector - visible bias parameters.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% Arranged [a, v], [a, v+1] ... [a+1, v] ...
% ---------------------------------


Params = zeros(NQSObj.Np,1);

W_shift = NQSObj.Wm.';
Params = Params + [NQSObj.av; NQSObj.bv; W_shift(:)]; % Will throw error if not same size.
end