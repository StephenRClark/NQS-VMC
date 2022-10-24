% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQS(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.Alpha = number of unique coupling sets or "hidden unit density".
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.av = (Nsl x 1) vector - visible bias parameters.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% Arranged [a, v], [a, v+1] ... [a+1, v] ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.
Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
for s = 1:Nsl
    if sum(NQSObj.OptInds(s,:)) ~= 0
        dLogp(s) = sum(Cfg_vec(SLInds==s)); % Insert d/da.
    end
end
dTheta = tanh(NQSObj.Theta);
% Accounting for shift structure of W matrix requires either construction
% of shifted Theta matrix or shifted Cfg vector - the latter is done here
for al=1:Alpha % Derivatives need to be computed by Alpha sector
    bInd = Nsl+al;
    if sum(NQSObj.OptInds(Nsl+al,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
    end
    for v = 1:Nv
        PInd = Nsl + Alpha + v + (al-1)*Nv;
        % For each layer labelled by a, find the indices of the associated translates.
        if sum(NQSObj.OptInds(PInd,:)) ~= 0
            for bd = 1:Ntr
                TInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                if VInd ~= 0
                    dLogp(PInd) = dLogp(PInd) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dW.
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
end