% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSOHTI(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim + Alpha + (VDim*Nv * Alpha).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% Properties added with translation invariance:
% - NQS.ati = (VDim x 1) vector - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x VDim*Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (VDim x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv*VDim x 1) for d/dW.
% Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Alpha = NQSObj.Alpha;
Ntr = numel(BondMap);
dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
OHVec = zeros(VDim,Nv);
for v = 1:VDim
    OHVec(v,:) = (Cfg_vec.' == NQSObj.VList(v));
    if NQSObj.OptInds(v) ~= 0
        dLogp(v) = sum(Cfg_vec==NQSObj.VList(v));
    end
end
OHVec = reshape(OHVec,VDim*Nv,1);
dTheta = tanh(NQSObj.Theta);
for a = 1:Alpha
    if NQSObj.OptInds(VDim+a) ~= 0
        dLogp(VDim+a) = sum(dTheta((1:Ntr)+(a-1)*Ntr));
    end
    for n = 1:Nv
        for v = 1:VDim
            PInd = VDim + Alpha + VDim*((n-1)+(a-1)*Nv) + v;
            if NQSObj.OptInds(PInd) ~= 0
                for b = 1:numel(BondMap)
                    HInd = b + (a-1)*Ntr; VInd = v + VDim*(BondMap{b}(n)-1);
                    dLogp(PInd) = dLogp(PInd) + OHVec(VInd)*dTheta(HInd);
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end