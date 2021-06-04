% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSMHAI(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.WH = (Nh x Nv) matrix - hidden-visible HH coupling terms.
% - NQS.WM = (Nh x Nv) matrix - hidden-visible MM coupling terms.
% - NQS.XH = (Nh x Nv) matrix - hidden-visible HM coupling terms.
% - NQS.XM = (Nh x Nv) matrix - hidden-visible MH coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% Properties added with translation invariance:
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.BHti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.BMti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.WHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.WMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
NQSObj.NsqVec = Cfg_vec.^2; % Log the initial configuration's square values.
NQSObj.Hv = (Cfg_vec == 0); % Create the vector that represents the holon operator.
NQSObj.Mv = (Cfg_vec - 1) .* (Cfg_vec > 0); % Create the vector that represents the multiplon operator.
% The multiplon operator is defined as M = n + H - 1.
NQSObj.ThetaH = (NQSObj.WH * NQSObj.Hv) + (NQSObj.XH * NQSObj.Mv) + NQSObj.BH; % Hidden holon effective angles.
NQSObj.ThetaM = (NQSObj.WM * NQSObj.Mv) + (NQSObj.XM * NQSObj.Hv) + NQSObj.BM; % Hidden multiplon effective angles.