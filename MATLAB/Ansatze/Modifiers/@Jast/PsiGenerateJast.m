% --- Exact normalised Jastrow amplitude generating function ---

function [Psi] = PsiGenerateJast(JastObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for Jastrow Modifier object:
% - Jast.N = number of sites (defined on input).
% - Jast.Np = number of variational Jastrow parameters.
% - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
% - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
% - Jast.Tj = (N x 1) vector - used to track on-site contributions.
% - Jast.JsV = (N x N) matrix - contains variational parameter indices for each site.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Js = JastObj.Js;
Psi = zeros(size(Basis,1),1);
for p = 1:numel(Psi)
    JfMat = Js .* (Basis(p,:).' * Basis(p,:));
    Psi(p) = exp(-0.5*sum(JfMat(:)));
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi)); 
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end