% --- General CPS wave function preparation function ---

% - Original by X Fang, updated by M Pei.

function [CPSObj] = PrepPsiCPSA(CPSObj,Cfg)
% This function initialises the CPS ansatz structure intermediate
% information (effective angles Theta and local indices VisInds) given an
% initial configuration.
% A indicates algebraic parameterisation.
% ---------------------------------
% Format for CPSA Modifier object:
% - CPSA.Nv = number of "visible" spins.
% - CPSA.Nh = number of "hidden" spins.
% - CPSA.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPSA.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPSA.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPSA.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPSA.HDim = 3 - this version features fixed hidden unit dimension.
% - CPSA.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPSA.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPSA.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPSA.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPSA.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% ---------------------------------

Nv = CPSObj.Nv; Nh = CPSObj.Nh; HDim = CPSObj.HDim;
Ind0 = CPSObj.Ind0; IndV = CPSObj.IndV;
CPSObj.VisInds = CPSObj.FullCfg(Cfg) + Ind0;
CPSObj.Theta = CPSObj.b;
for h = 1:Nh
    for v = 1:Nv
        if IndV(CPSObj.VisInds(v)) ~= 0
            for hd = 1:(HDim-1)
                CPSObj.Theta(h,hd) = CPSObj.Theta(h,hd) * CPSObj.W(IndV(CPSObj.VisInds(v)),hd,v,h);
            end
        end
    end
end
