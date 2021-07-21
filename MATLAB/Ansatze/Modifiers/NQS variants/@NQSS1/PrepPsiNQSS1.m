% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSS1(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
NQSObj.VisVec = Cfg_vec;
NQSObj.NsqVec = Cfg_vec.^2; % Log the initial configuration's square values. 
% Compute the effective angle for this configuration state.
NQSObj.Theta = NQSObj.b + (NQSObj.w*Cfg_vec) + (NQSObj.W*NQSObj.NsqVec);