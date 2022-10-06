% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSB(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis (which should be bosonic for a number-hidden
% NQS). This function will likely run into memory problems unless the
% number of sites is small.
% ---------------------------------
% Format for NQSB Modifier:
% - NQSB.Nv = number of "visible" units.
% - NQSB.Nh = number of "hidden" units.
% - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSB.a = (Nv x 1) vector - visible site bias.
% - NQSB.av = (Nsl x 1) vector - visible bias parameters.
% - NQSB.A = (Nv x 1) vector - visible site square bias.
% - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSB.b = (Nh x 1) vector - hidden site bias.
% - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSB.B = (Nh x 1) vector- hidden site square bias.
% - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSB.HDim = dimension of the hidden units.
% - NQSB.Theta = (Nh x 1) vector - effective angles.
% - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
N_cfgs = size(Basis,1);
A = (NQSObj.A.').*ones(N_cfgs,1); a = (NQSObj.a.').*ones(N_cfgs,1); 
Psi = exp(sum(A.*(Basis.^2) + (a.*Basis),2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum(Basis .* NQSObj.W(h,:),2) + NQSObj.b(h);
    B = NQSObj.B(h).*ones(N_cfgs,1);
    Psi = Psi .* NHTrace(Theta,B,NQSObj.HDim);
end
if ~isinf(max(abs(Psi)))
    Psi = Psi/max(abs(Psi)); % Pre-normalisation to avoid runaway arguments.
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end