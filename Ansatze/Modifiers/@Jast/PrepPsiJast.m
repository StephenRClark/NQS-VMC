% --- General Jastrow wave function preparation function ---

function [JastObj] = PrepPsiJast(JastObj,Cfg)
% This function initialises the Jastrow Modifier object with intermediate
% information (local Jastrow contributions) given an initial configuration.
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

% Construct vector representation of Cfg.
Cfg_vec = JastObj.FullCfg(Cfg);
% Sum over 2nd index accounts for fermionic representations.

% Calculate the local contribution to the overall Jastrow factor for each
% site in Cfg_vec. Can be done at once with a matrix multiplication.
JastObj.Tj = JastObj.Js * Cfg_vec;