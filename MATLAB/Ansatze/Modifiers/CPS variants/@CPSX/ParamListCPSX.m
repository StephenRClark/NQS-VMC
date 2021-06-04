% --- General NQS wave function parameter listing function ---

% - Original by X Fang, updated by M Pei.

function [Params] = ParamListCPSX(CPSObj)
% This function lists all the associated parameters of a general CPS in a
% single vector.
% X indicates exponential parameterisation.
% ---------------------------------
% Format for CPSX Modifier object:
% - CPSX.Nv = number of "visible" spins.
% - CPSX.Nh = number of "hidden" spins.
% - CPSX.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPSX.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPSX.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPSX.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPSX.HDim = 3 - this version features fixed hidden unit dimension.
% - CPSX.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPSX.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPSX.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPSX.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPSX.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% ---------------------------------
Nv = CPSObj.Nv; Nh = CPSObj.Nh; VDim = CPSObj.VDim; HDim = CPSObj.HDim;

Params = zeros(CPSObj.Np,1);

Params(1:(Nv * (VDim-1))) = CPSObj.a(:);
Params((Nv * (VDim-1)) + (1:(Nh * (HDim-1)))) = CPSObj.b(:);
Params((Nv * (VDim-1)) + (Nh * (HDim-1)) + (1:(Nh * Nv * (HDim-1) * (VDim-1)))) = CPSObj.W(:);
end