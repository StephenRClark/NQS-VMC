% --- General fermionic Pfaffian wave function update function ---

function PfafObj = PsiUpdatePfaf(PfafObj,P)
% This function updates the variational parameters of the fermionic
% Pfaffian ansatz from a vector of parameters P.
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
% Format for dLogp vector is a (Np x 1) vector of parameter derivatives.
% ---------------------------------

PfafObj.PfVar = PfafObj.PfVar + P;

cap = PfafObj.ParamCap;

% Sanity check the values of the ansatz:
PfafObj.PfVar(isinf(PfafObj.PfVar)) = 0;
PfafObj.PfVar(isnan(PfafObj.PfVar)) = 0;
ind = abs(PfafObj.PfVar)>cap;
PfafObj.PfVar(ind) = sign(PfafObj.PfVar(ind))*cap;

% Repackage PfFull with new values - easily performed using PfV, which is fixed:
% Construct full pairing term matrix using PfVar and PfV.
PairMat = PfafObj.PfVar(PfafObj.PfV);
PfafObj.PairMat = PairMat;