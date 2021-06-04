% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSMH(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
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
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/dBH.
% - (Nh x 1) for d/dBM.
% - (Nh*Nv x 1) for d/dW.
% - (Nh*Nv x 1) for d/dX.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.
ThetaH = NQSObj.ThetaH; ThetaM = NQSObj.ThetaM; HDim = NQSObj.HDim;

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
Hv = (Cfg_vec == 0); Mv = (Cfg_vec - 1) .* (Cfg_vec > 0); % Build the holon / multiplon vectors.

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

dLogp(1:Nv) = Cfg_vec; % Insert d/da.
dLogp((1:Nv)+Nv) = Cfg_vec.^2; % Insert d/dA.

dLogp((1:Nh)+2*Nv) = dTH_MHTrace(ThetaH,ThetaM,HDim); % Insert d/dBH.
dLogp((1:Nh)+2*Nv+Nh) = dTM_MHTrace(ThetaH,ThetaM,HDim); % Insert d/dBM.

dLogp((1:(Nv*Nh))+2*(Nv+Nh)) = reshape( ((dTH_MHTrace(Theta,ThetaM,HDim)*Hv.') + ...
    (dTM_MHTrace(Theta,ThetaM,HDim)*Mv.')), Nh*Nv, 1); % Insert d/dW.

dLogp((1:(Nv*Nh))+2*(Nv+Nh)+(Nv*Nh)) = reshape( ((dTH_MHTrace(Theta,ThetaM,HDim)*Mv.') + ...
    (dTM_MHTrace(Theta,ThetaM,HDim)*Xv.')), Nh*Nv, 1);% Insert d/dX.

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;