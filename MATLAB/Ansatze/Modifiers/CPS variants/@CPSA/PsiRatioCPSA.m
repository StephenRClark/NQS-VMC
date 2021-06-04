% --- General NQS wave function amplitude ratio function ---

% - Original by X Fang, updated by M Pei.

function [Ratio,Update] = PsiRatioCPSA(CPSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
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
% Format for Update is a vector of new effective angles ThetaP and an
% updated local index vector VisInds.
% ---------------------------------

% Set up new index vector IndsP.
IndV = CPSObj.IndV; VisIndsP = CPSObj.VisInds; Ind0 = CPSObj.Ind0; HDim = CPSObj.HDim;
for d = 1:Diff.num
    VisIndsP(Diff.pos(d)) = CPSObj.VisInds(Diff.pos(d)) + Diff.val(d);
end
ThetaP = CPSObj.Theta;
for h = 1:CPSObj.Nh
    for k = 1:Diff.num
        % Apply changes to each Theta vector.
        for hd = 1:(HDim-1)
            if CPSObj.IndV(CPSObj.VisInds(Diff.pos(k))) ~= 0
                ThetaP(h,hd) = ThetaP(h,hd) / CPSObj.W(IndV(CPSObj.VisInds(Diff.pos(k))),hd,Diff.pos(k),h);
            end
            if CPSObj.IndV(VisIndsP(Diff.pos(k))) ~= 0
                ThetaP(h,hd) = ThetaP(h,hd) * CPSObj.W(IndV(VisIndsP(Diff.pos(k))),hd,Diff.pos(k),h);
            end
        end
    end
end
Ratio = 1;
% Apply visible correlator part of the ratio.
for k = 1:Diff.num
    if CPSObj.IndV(CPSObj.VisInds(Diff.pos(k))) ~= 0
        Ratio = Ratio / CPSObj.a(Diff.pos(k),IndV(CPSObj.VisInds(Diff.pos(k))));
    end
    if CPSObj.IndV(VisIndsP(Diff.pos(k))) ~= 0
        Ratio = Ratio * CPSObj.a(Diff.pos(k),IndV(VisIndsP(Diff.pos(k))));
    end
end
% Perform the trace over the hidden dimension and take the product of trace
% ratios.
Ratio = Ratio * prod((1+sum(ThetaP,2))./(1+sum(Theta,2)));
Update.Theta = ThetaP;
Update.VisInds = VisIndsP;
end