% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQS(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
a = (NQSObj.a.').*ones(size(Basis,1),1);
Psi = exp(sum(a.*Basis,2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum(Basis .* NQSObj.W(h,:),2) + NQSObj.b(h);
    Psi = Psi .* cosh(Theta);
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end