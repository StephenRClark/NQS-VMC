% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSC(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQSC Modifier:
% - NQSC.Nv = number of "visible" units.
% - NQSC.Nh = number of "hidden" units.
% - NQSC.Np = number of parameters in the ansatz = 2*Alpha + 2*Alpha*Nv + 2*Nsl.
% - NQSC.a = (Nv x 1) vector - visible site bias.
% - NQSC.av = (Nsl x 1) vector - visible bias parameters.
% - NQSC.A = (Nv x 1) vector - visible site square bias.
% - NQSC.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSC.b = (Nh x 1) vector - hidden site bias.
% - NQSC.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSC.B = (Nh x 1) vector- hidden site square bias.
% - NQSC.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSC.w = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSC.wm = (Alpha x Nv) matrix - coupling parameters
% - NQSC.W = (Nh x Nv) matrix - hidden-square-visible coupling terms.
% - NQSC.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSC.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSC.HDim = dimension of the hidden units.
% - NQSC.Theta = (Nh x 1) vector - effective angles.
% - NQSC.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - vector of new effective angles ThetaP.
% Update.NsqVec - vector of new squared visible occupancies.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a configuration vector.
NQSObj.Theta = NQSObj.b + NQSObj.w*Cfg_vec; % Compute the effective angle for this configuration state.
NQSObj.NsqVec = Cfg_vec.^2; NQSObj.Theta = NQSObj.Theta + NQSObj.W*NQSObj.NsqVec; % Log the initial configuration's square values.
end