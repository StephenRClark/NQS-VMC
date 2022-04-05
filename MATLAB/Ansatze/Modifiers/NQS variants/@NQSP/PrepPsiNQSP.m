% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSP(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
% interactions. Max value VDim-1.
% - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
% interactions. Max value HDim-1.
% - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
% HOrder) + (Nv x VOrder)(Nh x HOrder)
% - NQS.a = (Nv x VOrder) matrix - visible site biases.
% - NQS.b = (Nh x HOrder) matrix - hidden site bias.
% - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.VisVec = (Nv x 1) vector - visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - matrix of new effective angles ThetaP.
% Update.VisVec - vector of new visible occupancies.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg)/(NQSObj.VDim-1); % Convert Cfg structure into a rescaled configuration vector.
Cfg_pow = reshape(Cfg_vec,1,NQSObj.Nv).^(reshape((1:NQSObj.VOrder),1,1,1,4)); % 1 x Nv x 1 x VOrder 
NQSObj.Theta = NQSObj.b + reshape(sum(sum(NQSObj.W.*Cfg_pow,2),4),NQSObj.Nh,NQSObj.HOrder); % Compute the effective angle for this configuration state.
end