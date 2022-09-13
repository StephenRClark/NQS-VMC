% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSU(NQSObjNew,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" spins.
% - NQSU.Nh = number of "hidden" spins.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
Nv = NQSObjNew.Nv; Nmax = NQSObjNew.VDim-1;
Psi = ones(size(Basis,1),1); % Visible bias contributions handled here.
for b = 1:size(Basis,1)
    UVec = zeros(Nmax,Nv);
    for v = 1:Nmax
        UVec(v,:) = (Basis(b,:) == NQSObjNew.VList(v+1));
    end
    UVec = reshape(UVec,Nmax*Nv,1);
    Theta = NQSObjNew.b + NQSObjNew.W*UVec;
    Psi(b) = exp(sum(NQSObjNew.a.*UVec)) * prod(cosh(Theta(:)));
end
if isinf(max(abs(Psi))) == 0
    Psi = Psi/max(abs(Psi)); % Pre-normalisation to avoid runaway arguments:
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end