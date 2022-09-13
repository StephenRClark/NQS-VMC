% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSA(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSA Modifier:
% - NQSA.Nv = number of "visible" spins.
% - NQSA.Nh = number of "hidden" spins.
% - NQSA.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSA.a = (Nv x 1) vector - visible site bias.
% - NQSA.av = (Nsl x 1) vector - visible bias parameters.
% - NQSA.A = (Nv x 1) vector - visible site square bias.
% - NQSA.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSA.b = (Nh x 1) vector - hidden site bias.
% - NQSA.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSA.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSA.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSA.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSA.Theta = (Nh x 1) vector - effective angles.
% - NQSA.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
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

dTheta = tanh(NQSObj.Theta);
for a = 1:Alpha
    bInd = 2*Nsl + a;
    if sum(NQSObj.OptInds(bInd,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(a-1)*Ntr));
    end
    for v = 1:Nv
        PInd = 2*Nsl + Alpha + (a-1)*Nv + v;
        if sum(NQSObj.OptInds(PInd,:)) ~= 0
            for b = 1:numel(BondMap)
                HInd = b + (a-1)*Ntr; VInd = BondMap{b}(v);
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