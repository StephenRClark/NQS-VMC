% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSMHAI(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS multiplon-holon ansatz, for a
% configuration specifed by the structure Cfg.
% NB: Translation invariance here assumes Nh/Nv integer and is currently
% only supported for hypercubic graphs.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.WH = (Nh x Nv) matrix - hidden-visible HH coupling terms.
% - NQS.WM = (Nh x Nv) matrix - hidden-visible MM coupling terms.
% - NQS.XH = (Nh x Nv) matrix - hidden-visible HM coupling terms.
% - NQS.XM = (Nh x Nv) matrix - hidden-visible MH coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% Properties added with translation invariance:
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.BHti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.BMti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.WHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.WMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (1 x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB.
% - (Alpha*Nv x 1) for d/dWH.
% - (Alpha*Nv x 1) for d/dWM.
% - (Alpha*Nv x 1) for d/dXH.
% - (Alpha*Nv x 1) for d/dXM.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Alpha = NQSObj.Alpha; % Hidden unit density, needs to be integer.
HDim = NQSObj.HDim; % Hidden unit dimension.

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
ThetaH = NQSObj.ThetaH; ThetaM = NQSObj.ThetaM;

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

if NQSObj.OptInds(1) == 1
    dLogp(1) = sum(Cfg_vec); % Insert d/da.
end
if NQSObj.OptInds(2) == 1
    dLogp(2) = sum(Cfg_vec.^2); % Insert d/dA
end
% Accounting for shift structure of W/X matrix requires either construction
% of shifted Theta matrix or shifted Cfg vector - the latter is done here
for a=1:Alpha % Derivatives need to be computed by Alpha sector
    if NQSObj.OptInds(2+a) == 1
        dLogp(2+a) = sum(dTH_MHTrace(ThetaH((1:Ntr)+(a-1)*Ntr),...
            ThetaM((1:Ntr)+(a-1)*Ntr),HDim)); % Insert d/dBH.
    end
    if NQSObj.OptInds(2+Alpha+a) == 1
        dLogp(2+Alpha+a) = sum(dTM_MHTrace(ThetaH((1:Ntr)+(a-1)*Ntr),...
            ThetaM((1:Ntr)+(a-1)*Ntr),HDim)); % Insert d/dBH.
    end
    for v = 1:Nv
        PIndWH = 2 + 2*Alpha + v + (a-1)*Nv; PIndWM = PIndWH + Alpha*Nv;
        PIndXH = PIndWM + Alpha*Nv; PIndXM = PIndXH + Alpha*Nv;
        % For each layer labelled by a, find the indices of the associated translates.
        for b = 1:Ntr
            TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
            if VInd ~= 0
                if NQSObj.OptInds(PIndWH) == 1
                dLogp(PIndWH) = dLogp(PIndWH) + (NQSObj.Hv(VInd)*dTH_MHTrace(ThetaH(TInd),...
                    ThetaM(TInd),HDim)); % Insert d/dWH.
                end
                if NQSObj.OptInds(PIndWM) == 1
                dLogp(PIndWM) = dLogp(PIndWM) + (NQSObj.Mv(VInd)*dTM_MHTrace(ThetaH(TInd),...
                    ThetaM(TInd),HDim)); % Insert d/dWH.
                end
                if NQSObj.OptInds(PIndXH) == 1
                dLogp(PIndXH) = dLogp(PIndXH) + (NQSObj.Mv(VInd)*dTH_MHTrace(ThetaH(TInd),...
                    ThetaM(TInd),HDim)); % Insert d/dWH.
                end
                if NQSObj.OptInds(PIndXM) == 1
                dLogp(PIndXM) = dLogp(PIndXM) + (NQSObj.Hv(VInd)*dTM_MHTrace(ThetaH(TInd),...
                    ThetaM(TInd),HDim)); % Insert d/dWH.
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
