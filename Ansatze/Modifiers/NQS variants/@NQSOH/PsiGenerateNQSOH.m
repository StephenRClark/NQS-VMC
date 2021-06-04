% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSOH(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Nv = NQSObj.Nv; VDim = NQSObj.VDim;
Psi = ones(size(Basis,1),1); % Visible bias contributions handled here.
for b = 1:size(Basis,1)
    OHVec = zeros(VDim,Nv);
    for v = 1:VDim
        OHVec(v,:) = (Basis(b,:) == NQSObj.VList(v));
    end
    OHVec = reshape(OHVec,VDim*Nv,1);
    Theta = NQSObj.b + NQSObj.W*OHVec;
    Psi(b) = exp(sum(NQSObj.a.*OHVec)) * prod(cosh(Theta(:)));
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end