% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSS1(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/db.
% - (Nh*Nv x 1) for d/dw.
% - (Nh*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
Cfg_sqr = Cfg_vec.^2;

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

dTheta = tanh(NQSObj.Theta);

dLogp(1:Nv) = Cfg_vec; % Insert d/da.
dLogp((1:Nv)+Nv) = Cfg_sqr; % Insert d/dA.
dLogp((1:Nh)+2*Nv) = dTheta; % Insert d/db.
dLogp((1:(Nv*Nh))+2*Nv+Nh) = reshape((dTheta*Cfg_vec.'),Nh*Nv,1); % Insert d/dw.
dLogp((1:(Nv*Nh))+2*Nv+Nh+Nv*Nh) = reshape((dTheta*Cfg_sqr.'),Nh*Nv,1); % Insert d/dW.

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end