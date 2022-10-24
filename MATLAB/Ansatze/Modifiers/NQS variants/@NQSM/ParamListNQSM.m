% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSM(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (3*Nsl x 1) for d/da.
% Arranged [sl, H], [sl, D], [sl, M], [sl+1, H], ...
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% Arranged [a, v], [a, v+1] ... [a+1, v], ...
% - (Alpha*Nv x 1) for d/dX.
% Arranged [a, v], [a, v+1] ...
% ---------------------------------
Params = zeros(NQSObj.Np,1);

W_shift = NQSObj.Wm.'; X_shift = NQSObj.Xm.';
Params = Params + [NQSObj.av; NQSObj.bv; W_shift(:); X_shift(:)]; % Will throw error if not same length.
end