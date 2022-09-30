% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSA(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQSA Modifier:
% - NQSA.Nv = number of "visible" units.
% - NQSA.Nh = number of "hidden" units.
% - NQSA.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSA.a = (Nv x 1) vector - visible site bias.
% - NQSA.av = (Nsl x 1) vector - visible bias parameters.
% - NQSA.A = (Nv x 1) vector - visible site square bias.
% - NQSA.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSA.b = (Nh x 1) vector - hidden site bias.
% - NQSA.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSA.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSA.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSA.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSA.Theta = (Nh x 1) vector - effective angles.
% - NQSA.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

Params = zeros(NQSObj.Np,1);

W_shift = NQSObj.Wm.';
Params = Params + [NQSObj.av; NQSObj.Av; NQSObj.b; W_shift(:)]; % Will throw error if not same length.
end