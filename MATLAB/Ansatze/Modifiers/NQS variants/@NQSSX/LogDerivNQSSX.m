% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSSX(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/db.
% - (Nh x 1) for d/dB.
% - (Nh*Nv x 1) for d/dW.
% - (Nh*Nv x 1) for d/dX.
% ---------------------------------

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
Cfg_sqr = Cfg_vec.^2;
dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
dTheta = dL_SqTrace(NQSObj.Theta,NQSObj.ThetaSq,NQSObj.HVals);
dThetaSq = dS_SqTrace(NQSObj.Theta,NQSObj.ThetaSq,NQSObj.HVals);
dLogp(1:NQSObj.Nv) = Cfg_vec; % Insert d/da.
dLogp((1:NQSObj.Nv)+NQSObj.Nv) = Cfg_sqr; % Insert d/dA.
dLogp((1:NQSObj.Nh)+2*NQSObj.Nv) = dTheta; % Insert d/db.
dLogp((1:NQSObj.Nh)+2*NQSObj.Nv+NQSObj.Nh) = dThetaSq; % Insert d/db.
dLogp((1:(NQSObj.Nv*NQSObj.Nh))+2*NQSObj.Nv+2*NQSObj.Nh) = ...
    reshape((dTheta*Cfg_vec.'),NQSObj.Nh*NQSObj.Nv,1); % Insert d/dW.
dLogp((1:(NQSObj.Nv*NQSObj.Nh))+2*NQSObj.Nv+2*NQSObj.Nh+NQSObj.Nv*NQSObj.Nh) = ...
    reshape((dThetaSq*Cfg_sqr.'),NQSObj.Nh*NQSObj.Nv,1); % Insert d/dX.
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end