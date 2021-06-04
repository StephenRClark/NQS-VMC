% --- General Jastrow wave function update function ---

function [JastObj] = PsiCfgUpdateJast(JastObj,Update)
% This function updates the intermediate configuration state information
% (on-site Jastrow contributions) retained in the Jastrow ansatz.
% ---------------------------------
% Format for Jastrow Modifier object:
% - Jast.N = number of sites (defined on input).
% - Jast.Np = number of variational Jastrow parameters.
% - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
% - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
% - Jast.Tj = (N x 1) vector - used to track on-site contributions.
% - Jast.JsV = (N x N) matrix - contains variational parameter indices for each site.
% ---------------------------------
% N.B: Jastrow factors here are assumed symmetric i.e.
% Js(i,j) = JS(j,i).
% ---------------------------------
% Format for Update is a vector of new local Jastrow contributions.
% ---------------------------------

% Just overwrite the on-site Jastrow information computed earlier.
JastObj.Tj = Update;
