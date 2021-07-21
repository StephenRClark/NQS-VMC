% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSSX(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis (which should be bosonic for a number-hidden
% NQS). This function will likely run into memory problems unless the
% number of sites is small.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Ab = (NQSObj.A.').*ones(size(Basis,1),1); ab = (NQSObj.a.').*ones(size(Basis,1),1);
Psi = exp(sum(Ab.*(Basis.^2) + (ab.*Basis),2)); % Visible bias contributions handled here.
for h = 1:NQSObj.Nh
    Theta = sum((Basis .* NQSObj.W(h,:)),2) + NQSObj.b(h);
    ThetaSq = sum((Basis.^2 .* NQSObj.X(h,:)),2) + NQSObj.B(h);
    Psi = Psi .* SqTrace(Theta,ThetaSq,NQSObj.HVal);
end
ModPsi = sqrt(sum(abs(Psi).^2)); % Try normalising every hidden unit to avoid numerical overflow.
if ModPsi == 0
    ModPsi = 1;
end
Psi = Psi/ModPsi;
end