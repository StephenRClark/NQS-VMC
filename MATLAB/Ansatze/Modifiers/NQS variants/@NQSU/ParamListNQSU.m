% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSU(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" spins.
% - NQSU.Nh = number of "hidden" spins.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nmax*Nsl x 1) for d/da.
% Arranged [sl, vd], [sl, vd+1], ... , [sl+1, vd], ...
% - (Alpha x 1) for d/db.
% - (Alpha*Nv*Nmax x 1) for d/dW.
% Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
% ---------------------------------

Nv = NQSObj.Nv; Alpha = NQSObj.Alpha; Nmax = NQSObj.VDim-1; Nsl = max(NQSObj.Graph.SLInds);

Params = zeros(NQSObj.Np,1);

Params(1:(Nmax*Nsl)) = NQSObj.av;
Params((1:Alpha)+Nmax*Nsl) = NQSObj.bv;
Params((1:(Alpha*Nmax*Nv))+Alpha+Nmax*Nsl) = reshape(NQSObj.Wm.',Alpha*Nmax*Nv,1);
end