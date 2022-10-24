% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSS1(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQSS1.Nv = number of "visible" spins.
% - NQSS1.Nh = number of "hidden" spins.
% - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
% - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSS1.a = (Nv x 1) vector - visible site bias.
% - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
% - NQSS1.A = (Nv x 1) vector - visible site square bias.
% - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSS1.b = (Nh x 1) vector - hidden site bias.
% - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
% - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
% - NQSS1.Theta = (Nh x 1) vector - effective angles.
% - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
NQSObj.VisVec = Cfg_vec;
NQSObj.NsqVec = Cfg_vec.^2; % Log the initial configuration's square values. 
% Compute the effective angle for this configuration state.
NQSObj.Theta = NQSObj.b + (NQSObj.w*Cfg_vec) + (NQSObj.W*NQSObj.NsqVec);
end