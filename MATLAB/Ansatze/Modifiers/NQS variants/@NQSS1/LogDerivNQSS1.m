% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSS1(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQSS1.Nv = number of "visible" spins.
% - NQSS1.Nh = number of "hidden" spins.
% - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
% - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSS1.a = (Nv x 1) vector - visible site bias.
% - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
% - NQSS1.A = (Nv x 1) vector - visible site square bias.
% - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSS1.b = (Nh x 1) vector - hidden site bias.
% - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
% - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
% - NQSS1.Theta = (Nh x 1) vector - effective angles.
% - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dw.
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.
Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
Cfg_sqr = Cfg_vec.^2;
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
for al=1:Alpha % Derivatives need to be computed by Alpha sector
    bInd = 2*Nsl + al;
    if sum(NQSObj.OptInds(bInd,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
    end
    for v = 1:NQSObj.Nv
        PIndw = 2*Nsl + Alpha + v + (al-1)*Nv;
        PIndW = PIndw + Alpha*Nv;
        % For each layer labelled by a, find the indices of the associated translates.
        for bd = 1:Ntr
            TInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
            if VInd ~= 0
                if sum(NQSObj.OptInds(PIndw,:)) ~= 0
                    dLogp(PIndw) = dLogp(PIndw) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dw.
                end
                if sum(NQSObj.OptInds(PIndW,:)) ~= 0
                    dLogp(PIndW) = dLogp(PIndW) + (Cfg_sqr(VInd)*dTheta(TInd)); % Insert d/dW.
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
end