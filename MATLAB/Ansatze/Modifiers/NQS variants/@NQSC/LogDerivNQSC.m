% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSC(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSC Modifier:
% - NQSC.Nv = number of "visible" units.
% - NQSC.Nh = number of "hidden" units.
% - NQSC.Np = number of parameters in the ansatz = 2*Alpha + 2*Alpha*Nv + 2*Nsl.
% - NQSC.a = (Nv x 1) vector - visible site bias.
% - NQSC.av = (Nsl x 1) vector - visible bias parameters.
% - NQSC.A = (Nv x 1) vector - visible site square bias.
% - NQSC.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSC.b = (Nh x 1) vector - hidden site bias.
% - NQSC.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSC.B = (Nh x 1) vector- hidden site square bias.
% - NQSC.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSC.w = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSC.wm = (Alpha x Nv) matrix - coupling parameters
% - NQSC.W = (Nh x Nv) matrix - hidden-square-visible coupling terms.
% - NQSC.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSC.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSC.HDim = dimension of the hidden units.
% - NQSC.Theta = (Nh x 1) vector - effective angles.
% - NQSC.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB
% - (Alpha*Nv x 1) for d/dw.
% - (Alpha*Nv x 1) for d/dW.
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
    if sum(NQSObj.OptInds(s,:))~=0
        dLogp(s) = sum(Cfg_vec(SLInds==s)); % Insert d/da.
    end
    if sum(NQSObj.OptInds(s+Nsl,:))~=0
        dLogp(s+Nsl) = sum(Cfg_vec(SLInds==s).^2); % Insert d/dA.
    end
end

dTheta = dT_NTrace(NQSObj.Theta,NQSObj.B,NQSObj.HDim);
dB = dB_NHTrace(NQSObj.Theta,NQSObj.B,NQSObj.HDim);
for al = 1:Alpha
    bInd = 2*Nsl + al; BInd = bInd + Alpha;
    if sum(NQSObj.OptInds(bInd,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
    end
    if sum(NQSObj.OptInds(BInd,:)) ~= 0
        dLogp(BInd) = sum(dB((1:Ntr)+(al-1)*Ntr)); % Insert d/dB.
    end
    for v = 1:Nv
        PIndw = 2*Nsl + 2*Alpha + (al-1)*Nv + v;
        PIndW = PIndw + Alpha*Nv;
        if sum(NQSObj.OptInds(PIndw,:)) ~= 0
            for bd = 1:numel(BondMap)
                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                if VInd ~= 0
                    dLogp(PIndw) = dLogp(PIndw) + Cfg_vec(VInd)*dTheta(HInd);
                end
            end
        end
        if sum(NQSObj.OptInds(PIndW,:)) ~= 0
            for bd = 1:numel(BondMap)
                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                if VInd ~= 0
                    dLogp(PIndW) = dLogp(PIndW) + (Cfg_vec(VInd)^2)*dTheta(HInd);
                end
            end
        end
    end
end

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;

end