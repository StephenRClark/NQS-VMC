% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSM(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.VDim = dimensions of the visible units.
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (3*Nsl x 1) for d/da.
% Arranged [sl, H], [sl, D], [sl, M], [sl+1, H], ...
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% Arranged [a, v], [a, v+1] ... [a+1, v], ...
% - (Alpha*Nv x 1) for d/dX.
% Arranged [a, v], [a, v+1] ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.
% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.
Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
OMatP = [(Cfg_vec==0), (Cfg_vec==2), (Cfg_vec>2)];
dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
for s = 1:Nsl
    for v = 1:3
        PInd = v + (s-1)*3;
        if sum(NQSObj.OptInds(PInd,:))~=0
            dLogp((1:3) + (s-1)*3) = sum(OMatP(SLInds==s,:),1).'; % Insert d/da
        end
    end
end
dTheta = tanh(NQSObj.Theta);
for al = 1:Alpha
    bInd = 3*Nsl + al;
    if sum(NQSObj.OptInds(bInd,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr));
    end
    for n = 1:NQSObj.Nv
        PIndW = 3*Nsl + Alpha + (al-1)*Nv + n;
        PIndX = PIndW + NQSObj.Alpha*NQSObj.Nv;
        % Big if / elseif to avoid extraneous if checking or loops
        if (sum(NQSObj.OptInds(PIndW,:))~=0) && (sum(NQSObj.OptInds(PIndX))~=0)
            for bd = 1:numel(BondMap)
                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                dLogp(PIndW) = dLogp(PIndW) + OMatP(VInd,1)*dTheta(HInd);
                dLogp(PIndX) = dLogp(PIndX) + OMatP(VInd,2)*dTheta(HInd);
            end
        elseif (sum(NQSObj.OptInds(PIndW,:))~=0) && (sum(NQSObj.OptInds(PIndX))==0)
            for bd = 1:numel(BondMap)
                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                dLogp(PIndW) = dLogp(PIndW) + OMatP(VInd,1)*dTheta(HInd);
            end
        elseif (sum(NQSObj.OptInds(PIndW,:))==0) && (sum(NQSObj.OptInds(PIndX))~=0)
            for bd = 1:numel(BondMap)
                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                dLogp(PIndX) = dLogp(PIndX) + OMatP(VInd,2)*dTheta(HInd);
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
end