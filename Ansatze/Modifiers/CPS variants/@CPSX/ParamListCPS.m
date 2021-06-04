% --- General NQS wave function parameter listing function ---

% - Original by X Fang, updated by M Pei.

function [Params] = ParamListCPS(CPSObj)
% This function lists all the associated parameters of a general CPS in a
% single vector.
% ---------------------------------
% Format for CPSX Modifier object:
% - CPS.Nv = number of "visible" spins.
% - CPS.Nh = number of "hidden" spins.
% - CPS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPS.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPS.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPS.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPS.HDim = 3 - this version features fixed hidden unit dimension.
% - CPS.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPS.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPS.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPS.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPS.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% ---------------------------------
Nv = CPSObj.Nv; Nh = CPSObj.Nh; VDim = CPSObj.VDim; HDim = CPSObj.HDim;

Params = zeros(CPSObj.Np,1);

Params(1:(Nv * (VDim-1))) = CPSObj.a(:);
Params((Nv * (VDim-1)) + (1:(Nh * (HDim-1)))) = CPSObj.b(:);
Params((Nv * (VDim-1)) + (Nh * (HDim-1)) + (1:(Nh * Nv * (HDim-1) * (VDim-1)))) = CPSObj.W(:);
end