% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSP(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
% interactions. Max value VDim-1.
% - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
% interactions. Max value HDim-1.
% - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
% HOrder) + (Nv x VOrder)(Nh x HOrder)
% - NQS.a = (Nv x VOrder) matrix - visible site biases.
% - NQS.b = (Nh x HOrder) matrix - hidden site bias.
% - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.VisVec = (Nv x 1) vector - visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x VOrder) x 1 for d/da.
% - (Nh x HOrder) x 1 for d/db.
% - (Nh x Nv) x (HOrder x VOrder) for d/dW.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Basis = Basis / (NQSObj.VDim-1);
Psi = ones(size(Basis,1),1);
for vo = 1:NQSObj.VOrder
    avo = (NQSObj.a(:,vo).').*ones(size(Basis,1),1);
    Psi = exp(sum(avo.*(Basis.^vo),2)); % Visible bias contributions handled here.
end
for h = 1:NQSObj.Nh
    Theta = zeros(size(Basis,1),NQSObj.HOrder);
    HVals = (0:(NQSObj.HDim-1)/(NQSObj.HDim-1));
    for ho = 1:NQSObj.HOrder
        for vo = 1:NQSObj.VOrder
            Theta(:,ho) = sum((Basis.^vo) .* NQSObj.W(h,:,ho,vo),2) + NQSObj.b(h,ho);
        end
    end
    Psi = Psi .* ;
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end