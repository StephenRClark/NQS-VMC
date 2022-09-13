% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSNH(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + Nv*Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.B = (Nh x 1) vector - hidden site square bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/db.
% - (Nh x 1) for d/dB.
% - (Nh*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

dLogp(1:Nv) = Cfg_vec; % Insert d/da.
dLogp((1:Nv)+Nv) = Cfg_vec.^2; % Insert d/dA.
dLogp((1:(Nv*Nh))+2*(Nv+Nh)) = reshape((dTSBTrace(NQSObj.Theta,NQSObj.HDim,NQSObj.B)*Cfg_vec.'),Nh*Nv,1); % Insert d/dW.
dLogp((1:Nh)+2*Nv) = dTSBTrace(NQSObj.Theta,NQSObj.HDim,NQSObj.B); % Insert d/db.
dLogp((1:Nh)+2*Nv+Nh) = dBSBTrace(NQSObj.Theta,NQSObj.HDim,NQSObj.B); % Insert d/dB.
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
