% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSM(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP, new
% visible occupancies VisVec and new species matrix OMatP.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
OMat = [(Cfg_vec==0), (Cfg_vec==2), (Cfg_vec>2)];
NQSObj.OMat = OMat; NQSObj.VisVec = Cfg_vec;
NQSObj.Theta = NQSObj.b + [NQSObj.W, NQSObj.X]*([OMat(:,1);OMat(:,2)]); % Compute the effective angle for this configuration state.