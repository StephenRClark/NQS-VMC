% --- Exact normalised Jastrow amplitude generating function ---

function [Psi] = PsiGenerateNNMB(NNMBObj,Basis)
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
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
GMB = NNMBObj.GMB; NNList = NNMBObj.NNList;
Psi = zeros(size(Basis,1),1);
for p = 1:numel(Psi)
    Db0 = Basis(p,:)>=2; Hl0 = Basis(p,:)==0;
    Db = Db0; Hl = Hl0;
    for n = 1:size(NNList,2)
        DbNN = ones(1,size(NNList,1)); HlNN = ones(1,size(NNList,1));
        DbNN(NNList(:,n)>0) = 1-Db0(NNList(NNList(:,n)>0,n)); 
        HlNN(NNList(:,n)>0) = 1-Hl0(NNList(NNList(:,n)>0,n));
        Db = Db .* HlNN; Hl = Hl .* DbNN;
    end
    Psi(p) = exp(GMB*sum(Db + Hl));
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi));
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end