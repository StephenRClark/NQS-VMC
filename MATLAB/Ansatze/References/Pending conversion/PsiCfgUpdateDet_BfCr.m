% --- General fermionic determinant wave function update function ---

function [Ansatz] = PsiCfgUpdateDet_BfCr(Update,Ansatz)
% This function updates the intermediate configuration state information
% (matrices) retained in the fermionic determinant ansatz addition.
% BfCr - backflow correlations, will alter determinant according to
% doublon-holon positions.
% ---------------------------------
% Format for fermionic determinant wavefunction addition:
% - Ansatz.Nf = number of fermions - set to Nv for spin models, specified at initialisation.
% - Ansatz.UFull = (2Nv x 2Nv) matrix - contains all available single particle orbitals.
% - Ansatz.UFe = (2Nv x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
% - Ansatz.WFe = (2Nv x Nf) matrix - elements are used for determinants in PsiRatio.
% - Ansatz.FermLoc = (2Nv x 1) vector - details locations of fermions by index for sign tracking purposes.
% Backflow correlation alterations:
% - Ansatz.UInv (Nf x 2Nv) replaces WFe as it is more convenient to work with.
% ---------------------------------

% Just overwrite the information computed earlier.
Ansatz.FermLoc = Update.FermLoc;
Ansatz.UInv = Ansatz.UInv + Update.UInv;

