% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSS1(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis (which should be bosonic for a number-hidden
% NQS). This function will likely run into memory problems unless the
% number of sites is small.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Ab = (NQSObj.A.').*ones(size(Basis,1),1); ab = (NQSObj.a.').*ones(size(Basis,1),1);
Psi = exp(sum(Ab.*(Basis.^2) + (ab.*Basis),2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum((Basis .* NQSObj.w(h,:)) + (Basis.^2 .* NQSObj.W(h,:)),2) + NQSObj.b(h);
    Psi = Psi .* cosh(Theta);    
end
ModPsi = sqrt(sum(abs(Psi).^2)); % Try normalising every hidden unit to avoid numerical overflow.
    if ModPsi == 0
        ModPsi = 1;
    end
    Psi = Psi/ModPsi;
end