% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSP(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQSP Modifier:
% - NQSP.Nv = number of "visible" spins.
% - NQSP.Nh = number of "hidden" spins.
% - NQSP.Np = number of parameters in the ansatz = (Nsl x VOrder) + (Alpha x
% HOrder) + (Nv x VOrder)(Alpha x HOrder)
% - NQSP.VDim = dimension of the visible units.
% - NQSP.HDim = dimension of the hidden units.
% - NQSP.VOrder = highest power of visible unit interactions. Max value VDim-1.
% - NQSP.HOrder = highest power of hidden unit interactions. Max value HDim-1.
% - NQSP.a = (Nv x VOrder) matrix - visible site biases.
% - NQSP.av = (Nsl x VOrder) matrix - visible bias parameters
% - NQSP.b = (Nh x HOrder) matrix - hidden site bias.
% - NQSP.bv = (Alpha x HOrder) matrix - hidden bias parameters.
% - NQSP.W = (Nh x Nv x HOrder x VOrder) array - hidden-visible coupling terms.
% - NQSP.Wm = (Alpha x Nv x HOrder x VOrder) array - hidden-visible coupling parameters
% - NQSP.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSP.Theta = (Nh x HOrder) matrix - effective angles by hidden order.
% - NQSP.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSP.Rescale = flag for visible unit rescaling to [0 1] interval.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - matrix of new effective angles ThetaP.
% Update.VisVec - vector of new visible occupancies.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg)*((NQSObj.VDim-1)^(-NQSObj.Rescale)); % Convert Cfg structure into a rescaled configuration vector.
Cfg_pow = reshape(Cfg_vec,1,NQSObj.Nv).^(reshape((1:NQSObj.VOrder),1,1,1,4)); % 1 x Nv x 1 x VOrder
NQSObj.Theta = NQSObj.b + reshape(sum(sum(NQSObj.W.*Cfg_pow,2),4),NQSObj.Nh,NQSObj.HOrder); % Compute the effective angle for this configuration state.
end