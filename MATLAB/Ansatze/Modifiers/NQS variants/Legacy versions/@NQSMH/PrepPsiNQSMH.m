% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSMH(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nv + 2*Nh + 2*Nv.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
% - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
NQSObj.NsqVec = Cfg_vec.^2; % Log the initial configuration's square values.
NQSObj.Hv = (Cfg_vec == 0); % Create the vector that represents the holon operator.
NQSObj.Mv = (Cfg_vec - 1) .* (Cfg_vec > 0); % Create the vector that represents the multiplon operator.
% The multiplon operator is defined as M = n + H - 1.
NQSObj.ThetaH = (NQSObj.W * NQSObj.Hv) + (NQSObj.X * NQSObj.Mv) + NQSObj.BH; % Hidden holon effective angles.
NQSObj.ThetaM = (NQSObj.W * NQSObj.Mv) + (NQSObj.X * NQSObj.Hv) + NQSObj.BM; % Hidden multiplon effective angles.