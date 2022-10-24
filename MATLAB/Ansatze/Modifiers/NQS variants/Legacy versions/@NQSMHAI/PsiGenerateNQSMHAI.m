% --- Exact normalised NQS amplitude generating function ---

function [Psi] = PsiGenerateNQSMHAI(NQSObj,Basis)
% This function computes the full normalised wavefunction Psi for a
% supplied many-body Basis. This function will likely run into memory
% problems unless the number of sites is small.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.WH = (Nh x Nv) matrix - hidden-visible HH coupling terms.
% - NQS.WM = (Nh x Nv) matrix - hidden-visible MM coupling terms.
% - NQS.XH = (Nh x Nv) matrix - hidden-visible HM coupling terms.
% - NQS.XM = (Nh x Nv) matrix - hidden-visible MH coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% Properties added with translation invariance:
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.BHti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.BMti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.WHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.WMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------

% Basis should be a N_cfg x N matrix. Ensure visible biases align with
% configurations.
A = (NQSObj.A.').*ones(size(Basis,1),1); a = (NQSObj.a.').*ones(size(Basis,1),1);
Psi = exp(sum(A.*(Basis.^2) + (a.*Basis),2)); % Visible bias contributions handled here.
H = (Basis == 0); M = (Basis - 1) .* (Basis > 0);
for h = 1:NQSObj.Nh
    ThetaH = sum((H .* NQSObj.WH(h,:))+(M .* NQSObj.XH(h,:)),2) + NQSObj.BH(h);
    ThetaM = sum((M .* NQSObj.WM(h,:))+(H .* NQSObj.XM(h,:)),2) + NQSObj.BM(h);
    Psi = Psi .* MHTrace(ThetaH,ThetaM,NQSObj.HDim);
end
ModPsi = sqrt(sum(abs(Psi).^2));
Psi = Psi/ModPsi;
end