% --- General fermionic determinant wave function update function ---

function [SDetObj] = PsiCfgUpdateSDet(SDetObj,Update)
% This function updates the intermediate configuration state information
% (matrices) retained in the fermionic determinant ansatz addition.
% ---------------------------------
% Format for SDet Reference:
% - SDet.Nf = (1 x 2) vector - number of up/down fermions respectively.
% - SDet.Orbitals = (2N x 2N) matrix - contains all available single particle orbitals.
% - SDet.OrbMat = (2N x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
% - SDet.DetMat = (2N x Nf) matrix - elements are used for determinants in PsiRatio.
% - SDet.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
% - SDet.Np = number of variational parameters associated with SDet Reference.
% SDet properties used in variational version:
% - SDet.CArr = (2N x 2N x Np) array - connectivity array for the reference Hamiltonian.
% - SDet.WArr = (2N x 2N x Np) array - transformed connectivity array for the reference Hamiltonian.
% - SDet.EnFac = (2N x 2N) matrix - elements are used in LogDeriv function.
% - SDet.HVar = (Np x 1) vector - variational parameters in the reference Hamiltonian used.
% ---------------------------------
% Format for Update is a struct containing updates for DetMat and FermLoc.
% ---------------------------------

% Just overwrite the information computed earlier.
SDetObj.FermLoc = Update.FermLoc; N = numel(Update.FermLoc)/2;
SDetObj.Nf = [sum(SDetObj.FermLoc(1:N)>0), sum(SDetObj.FermLoc((1:N)+N)>0)];
SDetObj.DetMat = SDetObj.DetMat + Update.DetMatD;

