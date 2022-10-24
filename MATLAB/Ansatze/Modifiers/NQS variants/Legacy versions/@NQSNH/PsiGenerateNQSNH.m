% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSNH(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis (which should be bosonic for a number-hidden
% NQS). This function will likely run into memory problems unless the
% number of sites is small.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + Nv*Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.B = (Nh x 1) vector - hidden site square bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
A = (NQSObj.A.').*ones(size(Basis,1),1); a = (NQSObj.a.').*ones(size(Basis,1),1); 
Psi = exp(sum(A.*(Basis.^2) + (a.*Basis),2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum(Basis .* NQSObj.W(h,:),2) + NQSObj.b(h);
    B = NQSObj.B(h).*ones(size(Basis,1),1);
    Psi = Psi .* NHTrace(Theta,B,NQSObj.HDim);
end
if isinf(max(abs(Psi))) == 0
    Psi = Psi/max(abs(Psi));
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end