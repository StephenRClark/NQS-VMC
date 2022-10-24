% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSP(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQSP Modifier:
% - NQSP.Nv = number of "visible" units.
% - NQSP.Nh = number of "hidden" units.
% - NQSP.Np = number of parameters in the ansatz = (Nsl x VOrder) + (Alpha x
% HOrder) + (Nv x VOrder)(Alpha x HOrder)
% - NQSP.VDim = dimension of the visible units.
% - NQSP.HDim = dimension of the hidden units.
% - NQSP.VOrder = highest power of visible unit interactions. Max value VDim-1.
% - NQSP.HOrder = highest power of hidden unit interactions. Max value HDim-1.
% - NQSP.a = (Nv x VOrder) matrix - visible site biases.
% - NQSP.av = (Nsl x VOrder) matrix - visible bias parameters
% - NQSP.b = (Nh x HOrder) matrix - hidden site bias.
% - NQSP.bv = (Alpha x HOrder) matrix - hidden bias parameters.
% - NQSP.W = (Nh x Nv x HOrder x VOrder) array - hidden-visible coupling terms.
% - NQSP.Wm = (Alpha x Nv x HOrder x VOrder) array - hidden-visible coupling parameters
% - NQSP.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSP.Theta = (Nh x HOrder) matrix - effective angles by hidden order.
% - NQSP.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSP.Rescale = flag for visible unit rescaling to [0 1] interval.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
N_cfgs = size(Basis,1);
Nmax = NQSObj.VDim - 1; VOrder = NQSObj.VOrder;
HDim = NQSObj.HDim; HOrder = NQSObj.HOrder; Rescale = NQSObj.Rescale;
Basis = Basis * (Nmax^(-Rescale));
Psi = ones(N_cfgs,1);
for vo = 1:VOrder
    avo = (NQSObj.a(:,vo).').*ones(N_cfgs,1);
    Psi = exp(sum(avo.*(Basis.^vo),2)); % Visible bias contributions handled here.
end
HVals = (0:(HDim-1)) * ((HDim-1)^(-Rescale));
for h = 1:NQSObj.Nh
    Theta = zeros(N_cfgs,NQSObj.HOrder); % Clear it out for each hidden unit
    ThetaPow = zeros(N_cfgs,HDim); % = sum_ho Theta_ho * h^ho 
    for ho = 1:HOrder
        for vo = 1:VOrder
            Theta(:,ho) = Theta(:,ho) + sum((Basis.^vo) .* NQSObj.W(h,:,ho,vo),2) + NQSObj.b(h,ho);
        end
        ThetaPow = ThetaPow + (Theta(:,ho).*(HVals.^ho)); % = Theta * h^ho
    end
    Psi = Psi .* sum(exp(ThetaPow),2);
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi)); % Pre-normalisation to avoid runaway arguments:
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end