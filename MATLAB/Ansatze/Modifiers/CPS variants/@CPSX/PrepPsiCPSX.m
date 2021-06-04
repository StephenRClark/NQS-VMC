% --- General CPS wave function preparation function ---

% - Original by X Fang, updated by M Pei.

function [CPSObj] = PrepPsiCPSX(CPSObj,Cfg)
% This function initialises the CPS ansatz structure intermediate
% information (effective angles Theta and local indices VisInds) given an
% initial configuration.
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

Nv = CPSObj.Nv; Nh = CPSObj.Nh; HDim = CPSObj.HDim;
Ind0 = CPSObj.Ind0; IndV = CPSObj.IndV;
CPSObj.VisInds = CPSObj.FullCfg(Cfg) + Ind0;
CPSObj.Theta = CPSObj.b;
for h = 1:Nh
    for v = 1:Nv
        if IndV(CPSObj.VisInds(v)) ~= 0
            for hd = 1:(HDim-1)
                CPSObj.Theta(h,hd) = CPSObj.Theta(h,hd) + CPSObj.W(IndV(CPSObj.VisInds(v)),hd,v,h);
            end
        end
    end
end
