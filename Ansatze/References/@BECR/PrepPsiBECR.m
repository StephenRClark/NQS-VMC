% --- General bosonic reference wave function preparation function ---

function [BECRObj] = PrepPsiBECR(BECRObj,Cfg)
% This function initialises the bosonic reference ansatz structure 
% intermediate information (occupations) given an initial configuration.
% ---------------------------------
% Format for Bose Einstein Condensate Reference:
% - BECR.Nb = number of bosons - assumed fixed for the most part.
% - BECR.SPO = (Nv x 1) vector - single particle boson orbital amplitudes.
% - BECR.Occ = (Nv x 1) vector - onsite boson occupation numbers.
% ---------------------------------

% Simply add configuration occupation numbers to Ansatz.
BECRObj.Occ = Cfg.occ;