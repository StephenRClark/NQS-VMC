% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSP(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
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
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x VOrder) x 1 for d/da.
% - (Nh x HOrder) x 1 for d/db.
% - (Nh x Nv) x (HOrder x VOrder) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.
HOrder = NQSObj.HOrder; VOrder = NQSObj.VOrder;

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
NQSObj.VisVec = Cfg_vec; % Ensure VisVec is properly assigned.

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

VisPow = Cfg_vec(:) .^ (1:VOrder);
HidPow = ((0:(NQSObj.HDim-1))/(NQSObj.HDim-1).') .^ (1:NQSObj.HOrder);
HidArray = reshape(HidPow.',1,HOrder,NQSObj.HDim);
dTheta = (sum(exp(NQSObj.Theta.*HidArray).*HidArray,3))./sum(exp(NQSObj.Theta.*HidArray),3);
for vo = 1:VOrder
    dLogp(1:Nv+(vo-1)*Nv) = VisPow(:,vo); % Insert d/da.
end
for ho = 1:HOrder
    dLogp(1:Nh+(ho-1)*Nh+Nv*VOrder) = dTheta(:,ho);
    for v = 1:Nv
        for vo = 1:VOrder
            IndP = Nv*VOrder + Nh*HOrder + Nh*(v-1) + Nh*Nv*(ho-1) + Nh*Nv*HOrder*(vo-1);
            dLogp(1:Nh+IndP) = dTheta(:,ho)*VisPow(v,vo); % Insert d/dW.
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = dLogp * (NQSObj.OptInds~=0);
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end