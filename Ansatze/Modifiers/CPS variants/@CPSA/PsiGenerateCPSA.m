% --- Exact normalised CPS amplitude generating function ---

% - Original by X Fang, updated by M Pei.

function [Psi] = PsiGenerateCPSA(CPSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
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

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.

Psi = ones(size(Basis,1),1);
Basis = Basis + CPSObj.Ind0;
for b = 1:size(Basis,1)
    VInds = CPSObj.IndV(Basis(b,:)); N0Inds = (VInds~=0);
    Theta = CPSObj.b;
    for v = 1:CPSObj.Nv % Have to assemble Theta by visible units.
        if N0Inds(v)
            Psi(b) = Psi(b) * CPSObj.a(v,VInds(v));
            Theta = Theta .* (reshape(CPSObj.W(VInds(v),:,v,:),(CPSObj.HDim-1),CPSObj.Nh).');
        end
    end
    Psi(b) = Psi(b) * prod((1+sum(Theta,2))./CPSObj.HDim);
end

ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end