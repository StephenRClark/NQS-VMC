% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQS(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.Alpha = number of unique coupling sets or "hidden unit density".
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.av = (Nsl x 1) vector - visible bias parameters.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a configuration vector.
NQSObj.Theta = NQSObj.b + NQSObj.W*Cfg_vec; % Compute the effective angle for this configuration state.
end