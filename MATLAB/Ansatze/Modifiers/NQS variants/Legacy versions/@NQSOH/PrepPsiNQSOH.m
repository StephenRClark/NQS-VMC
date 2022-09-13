% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSOH(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP and
% new one-hot vector OHVecP.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
OHVec = zeros(NQSObj.VDim,NQSObj.Nv);
for v = 1:NQSObj.VDim
    OHVec(v,:) = (Cfg_vec.' == NQSObj.VList(v));
end
OHVec = reshape(OHVec,NQSObj.VDim*NQSObj.Nv,1);
NQSObj.OHVec = OHVec;
NQSObj.Theta = NQSObj.b + NQSObj.W*OHVec; % Compute the effective angle for this configuration state.