% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSS1TI(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS Spin ansatz, for a configuration
% specifed by the structure Cfg.
% NB: Translation invariance here assumes Nh/Nv integer and is currently
% only supported for hypercubic graphs.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Alpha + 2 + Alpha.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% Properties added with translation invariance:
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (1 x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dw.
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
Cfg_sqr = Cfg_vec.^2;
dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
dTheta = tanh(NQSObj.Theta);
if NQSObj.OptInds(1) == 1
    dLogp(1) = sum(Cfg_vec); % Insert d/da.
end
if NQSObj.OptInds(2) == 1
    dLogp(2) = sum(Cfg_sqr); % Insert d/dA.
end
for a=1:NQSObj.Alpha % Derivatives need to be computed by Alpha sector
    if NQSObj.OptInds(2+a) == 1
        dLogp(2+a) = sum(dTheta((1:Ntr)+(a-1)*Ntr)); % Insert d/db.
    end
    for v = 1:NQSObj.Nv
        PIndw = 2 + NQSObj.Alpha + v + (a-1)*NQSObj.Nv;
        PIndW = PIndw + NQSObj.Alpha*NQSObj.Nv;
        % For each layer labelled by a, find the indices of the associated translates.
        for b = 1:Ntr
            TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
            if VInd ~= 0
                if NQSObj.OptInds(PIndw) == 1
                    dLogp(PIndw) = dLogp(PIndw) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dw.
                end
                if NQSObj.OptInds(PIndW) == 1
                    dLogp(PIndW) = dLogp(PIndW) + (Cfg_sqr(VInd)*dTheta(TInd)); % Insert d/dW.
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end