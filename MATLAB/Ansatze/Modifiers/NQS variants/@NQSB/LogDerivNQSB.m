% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSB(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSB Modifier:
% - NQSB.Nv = number of "visible" spins.
% - NQSB.Nh = number of "hidden" spins.
% - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSB.a = (Nv x 1) vector - visible site bias.
% - NQSB.av = (Nsl x 1) vector - visible bias parameters.
% - NQSB.A = (Nv x 1) vector - visible site square bias.
% - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSB.b = (Nh x 1) vector - hidden site bias.
% - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSB.B = (Nh x 1) vector- hidden site square bias.
% - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSB.HDim = dimension of the hidden units.
% - NQSB.Theta = (Nh x 1) vector - effective angles.
% - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;

Alpha = NQSObj.Alpha; % Number of unique coupling sets.
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

dTheta = dTSBTrace(NQSObj.Theta,NQSObj.HDim,NQSObj.B);
dB = dBSBTrace(NQSObj.Theta,NQSObj.HDim,NQSObj.B);
for al = 1:NQSObj.Alpha
    bInd = 2*Nsl + al; BInd = bInd + NQSObj.Alpha;
    if sum(NQSObj.OptInds(bInd,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
    end
    if sum(NQSObj.OptInds(BInd,:)) ~= 0
        dLogp(BInd) = sum(dB((1:Ntr)+(al-1)*Ntr)); % Insert d/dB.
    end
    for v = 1:NQSObj.Nv
        PInd = 2*Nsl + 2*NQSObj.Alpha + (al-1)*NQSObj.Nv + v;
        if sum(NQSObj.OptInds(PInd,:)) ~= 0
            for bd = 1:numel(BondMap)
                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                if VInd ~= 0
                    dLogp(PInd) = dLogp(PInd) + Cfg_vec(VInd)*dTheta(HInd);
                end
            end
        end
    end
end

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;

end