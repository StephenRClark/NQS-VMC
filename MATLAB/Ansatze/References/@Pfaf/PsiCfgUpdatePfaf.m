% --- General fermionic Pfaffian wave function update function ---

function [PfafObj] = PsiCfgUpdatePfaf(PfafObj,Update)
% This function updates the intermediate configuration state information
% (matrices) retained in the fermionic Pfaffian ansatz addition.
% ---------------------------------
% Format for fermionic Pfaffian Reference:
% - Pfaf.Nf = (1 x 2) vector - number of up/down fermions respectively.
% - Pfaf.PairMat = (2N x 2N) matrix - contains all pairing terms.
% - Pfaf.PfI = (Nf x Nf) matrix - inverse of reduced PfFull matrix.
% - Pfaf.PfG = (2N x Nf) matrix - matrix used for ratio calculations.
% - Pfaf.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
% - Pfaf.Np = number of variational parameters associated with Pfaf Reference.
% Pfaf properties used in variational version:
% - Pfaf.PfV = (2N x 2N) array - logs which variational parameters make up the elements of Pfaf.PairMat.
% - Pfaf.PfVR = (Nf x Nf) array - reduced matrix constructed from PfV.
% - Pfaf.PfVar = (Np x 1) vector - variational parameters in PfFull.
% ---------------------------------
% Format for Update is a struct containing updates for PfI, PfG and FermLoc.
% ---------------------------------

% Just overwrite the information computed earlier.
PfafObj.FermLoc = Update.FermLoc;
PfafObj.PfG = PfafObj.PfG + Update.PfG;
PfafObj.PfI = PfafObj.PfI + Update.PfI;

% LogDerivPfaf depends highly on PfV being correct - worth explicitly
% reconstructing, since it's only a variable call.
PfInds = zeros(sum(PfafObj.Nf),1);
for i = 1:sum(PfafObj.Nf)
    PfInds(i) = find(PfafObj.FermLoc == i);
end
PfafObj.PfVR = PfafObj.PfV(PfInds,PfInds);