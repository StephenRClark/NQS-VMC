% --- General CPS wave function random initialisation function ---

% - Original by X Fang, updated by M Pei.

function [CPSObj] = RandomInitPsiCPS(CPSObj,Params)
% This function populates random initial CPS ansatz structure. The input
% Ansatz is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for CPS Modifier object:
% - CPS.Nv = number of "visible" spins.
% - CPS.Nh = number of "hidden" spins.
% - CPS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPS.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPS.b = (Nh x (HDim-1_) matrix - hidden site vector elements.
% - CPS.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPS.HDim = 3 - this version features fixed hidden unit dimension.
% - CPS.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPS.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPS.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPS.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPS.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% ---------------------------------

Nv = CPSObj.Nv; VDim = CPSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = CPSObj.Nh; HDim = CPSObj.HDim; % Number of "hidden" spins and hidden dimension.

CPSObj.Np = (VDim-1)*Nv + (HDim-1)*Nh + ((VDim-1)*(HDim-1) * (Nv*Nh)); % The number of variational parameters.
CPSObj.a = zeros(Nv,VDim-1);
CPSObj.b = zeros(Nh,HDim-1);
CPSObj.W = zeros(VDim-1,HDim-1,Nv,Nh);
CPSObj.Theta = zeros(Nh,HDim-1);
for v = 1:(Nv * (VDim-1))
    CPSObj.a(v) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h=1:(Nh * (HDim-1))
    CPSObj.b(h) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for k = 1:(Nv * Nh * (VDim-1)*(HDim-1))
    CPSObj.W(k) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
CPSObj.OptInds = ones(CPSObj.Np,1);
