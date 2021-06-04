% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSOH(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
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
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (VDim*Nv x 1) for d/da.
% Arranged [v, vd], [v, vd+1], ... , [v+1, vd], ...
% - (Nh x 1) for d/db.
% - (Nh*Nv*VDim x 1) for d/dW.
% Arranged [h, v, vd], [h, v, vd+1], ... ,[h, v+1, vd], ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = NQSObj.Nh; % Number of "hidden" spins.

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
OHVec = zeros(VDim,Nv);
for v = 1:VDim
    OHVec(v,:) = (Cfg_vec.' == NQSObj.VList(v));
end
OHVec = reshape(OHVec,VDim*Nv,1);

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
dTheta = tanh(NQSObj.Theta);

dLogp(1:(VDim*Nv)) = OHVec; % Insert d/da.
dLogp((1:Nh)+(VDim*Nh)) = dTheta; % Insert d/db.
dLogp((1:(Nh*Nv*VDim))+Nh+(VDim*Nh)) = reshape((dTheta * OHVec.'),Nh*Nv*VDim,1);
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end