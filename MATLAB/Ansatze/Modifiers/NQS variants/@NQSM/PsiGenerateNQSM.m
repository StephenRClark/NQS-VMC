% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSM(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
N_cfgs = size(Basis,1);
Nv = NQSObj.Nv;
Psi = ones(N_cfgs,1); % Visible bias contributions handled here.
for b = 1:N_cfgs
    OVec = [(Basis(b,:).'==0);(Basis(b,:).'==2);(Basis(b,:).'>2)];
    Theta = NQSObj.b + [NQSObj.W, NQSObj.X]*OVec(1:(2*Nv));
    Psi(b) = exp(sum(NQSObj.a.*OVec)) * prod(cosh(Theta(:)));
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi)); % Pre-normalisation to avoid runaway arguments.
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end