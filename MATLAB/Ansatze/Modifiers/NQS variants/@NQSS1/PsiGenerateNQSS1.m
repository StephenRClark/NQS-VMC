% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSS1(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis (which should be bosonic for a number-hidden
% NQS). This function will likely run into memory problems unless the
% number of sites is small.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQSS1.Nv = number of "visible" spins.
% - NQSS1.Nh = number of "hidden" spins.
% - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
% - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSS1.a = (Nv x 1) vector - visible site bias.
% - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
% - NQSS1.A = (Nv x 1) vector - visible site square bias.
% - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSS1.b = (Nh x 1) vector - hidden site bias.
% - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
% - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
% - NQSS1.Theta = (Nh x 1) vector - effective angles.
% - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
N_cfgs = size(Basis,1);
Ab = (NQSObj.A.').*ones(N_cfgs,1); ab = (NQSObj.a.').*ones(N_cfgs,1);
Psi = exp(sum(Ab.*(Basis.^2) + (ab.*Basis),2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum((Basis .* NQSObj.w(h,:)) + (Basis.^2 .* NQSObj.W(h,:)),2) + NQSObj.b(h);
    Psi = Psi .* cosh(Theta);    
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi)); % Pre-normalisation to avoid runaway arguments.
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end