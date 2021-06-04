% --- Exact normalised BEC amplitude generating function ---

function [Psi] = PsiGenerateBECR(BECRObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for Bose Einstein Condensate Reference:
% - BECR.Nb = number of bosons - assumed fixed for the most part.
% - BECR.SPO = (Nv x 1) vector - single particle boson orbital amplitudes.
% - BECR.Occ = (Nv x 1) vector - onsite boson occupation numbers.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Psi = zeros(size(Basis,1),1); SPO = reshape(BECRObj.SPO,1,numel(BECRObj.SPO));
for p = 1:numel(Psi)
    Psi(p) = sqrt(factorial(sum(Basis(p,:)))/prod(factorial(Basis(p,:)))) ...
        * prod(SPO.^Basis(p,:));
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end