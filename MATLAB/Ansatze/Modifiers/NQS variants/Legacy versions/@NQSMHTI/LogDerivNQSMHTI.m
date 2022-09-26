% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSMHTI(NQSObj,Cfg)
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
% - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
% - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
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
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (1 x 1) for d/dA.
% - (Alpha x 1) for d/dBH.
% - (Alpha x 1) for d/dBM.
% - (Alpha*Nv x 1) for d/dW.
% - (Alpha*Nv x 1) for d/dX.
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

dThetaH = dTH_MHTrace(ThetaH,ThetaM,HDim);
dThetaM = dTM_MHTrace(ThetaH,ThetaM,HDim);

% Accounting for shift structure of W/X matrix requires either construction
% of shifted Theta matrix or shifted Cfg vector - the latter is done here
for a=1:Alpha % Derivatives need to be computed by Alpha sector
    if NQSObj.OptInds(2+a) == 1
        dLogp(2+a) = sum(dThetaH((1:Ntr)+(a-1)*Ntr)); % Insert d/dBH.
    end
    if NQSObj.OptInds(2+Alpha+a) == 1
        dLogp(2+Alpha+a) = sum(dThetaM((1:Ntr)+(a-1)*Ntr)); % Insert d/dBM.
    end
    for v = 1:Nv
        PIndW = 2 + 2*Alpha + v + (a-1)*Nv; PIndX = PIndW + Alpha*Nv;
        % For each layer labelled by a, find the indices of the associated translates.
        if NQSObj.OptInds(PIndW) == 1
            for b = 1:Ntr
                TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
                if VInd ~= 0
                    dLogp(PIndW) = dLogp(PIndW) + (NQSObj.Hv(VInd)*dThetaH(TInd)) ...
                        + (NQSObj.Mv(VInd)*dThetaM(TInd)); % Insert d/dW.
                end
            end
        end
        if NQSObj.OptInds(PIndX) == 1
            for b = 1:Ntr
                TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
                if VInd ~= 0
                    dLogp(PIndX) = dLogp(PIndX) + (NQSObj.Mv(VInd)*dThetaH(TInd)) ...
                        + (NQSObj.Hv(VInd)*dThetaM(TInd)); % Insert d/dX.
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
