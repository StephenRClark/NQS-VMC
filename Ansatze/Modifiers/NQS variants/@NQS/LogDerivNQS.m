% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQS(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nh x 1) for d/db.
% - (Nh*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

dLogp(1:Nv) = Cfg_vec; % Insert d/da.
dLogp((Nv+Nh+1):(Nv+Nh+Nv*Nh)) = reshape((tanh(NQSObj.Theta)*Cfg_vec.'),Nh*Nv,1); % Insert d/dW.
dLogp((Nv+1):(Nv+Nh)) = tanh(NQSObj.Theta); % Insert d/db.
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
