% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSA(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis (which should be bosonic for a number-hidden
% NQS). This function will likely run into memory problems unless the
% number of sites is small.
% ---------------------------------
% Format for NQSA Modifier:
% - NQSA.Nv = number of "visible" units.
% - NQSA.Nh = number of "hidden" units.
% - NQSA.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSA.a = (Nv x 1) vector - visible site bias.
% - NQSA.av = (Nsl x 1) vector - visible bias parameters.
% - NQSA.A = (Nv x 1) vector - visible site square bias.
% - NQSA.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSA.b = (Nh x 1) vector - hidden site bias.
% - NQSA.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSA.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSA.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSA.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSA.Theta = (Nh x 1) vector - effective angles.
% - NQSA.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
N_cfgs = size(Basis,1);
A = (NQSObj.A.').*ones(N_cfgs,1); a = (NQSObj.a.').*ones(N_cfgs,1); 
Psi = exp(sum(A.*(Basis.^2) + (a.*Basis),2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum(Basis .* NQSObj.W(h,:),2) + NQSObj.b(h);
    Psi = Psi .* cosh(Theta);
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi)); % Pre-normalisation to avoid runaway arguments:
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end