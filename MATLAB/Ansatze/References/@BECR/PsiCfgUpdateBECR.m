% --- General bosonic reference wave function update function ---

function [BECRObj] = PsiCfgUpdateBECR(BECRObj,OccP)
% This function updates the intermediate configuration state information
% (occupations) retained in the bosonic reference ansatz addition.
% ---------------------------------
% Format for Bose Einstein Condensate Reference:
% - BECR.Nb = number of bosons - assumed fixed for the most part.
% - BECR.SPO = (Nv x 1) vector - single particle boson orbital amplitudes.
% - BECR.Occ = (Nv x 1) vector - onsite boson occupation numbers.
% ---------------------------------
% Format for Update is a vector containing the new boson occupation numbers.
% ---------------------------------

% Just overwrite the information computed earlier.
BECRObj.Occ = OccP;

