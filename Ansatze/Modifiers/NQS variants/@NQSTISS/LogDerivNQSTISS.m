% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSTISS(NQSObj,HilbertObj,GraphObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS Spin ansatz, for a configuration
% specifed by the structure Cfg.
% NB: Translation invariance here assumes Nh/Nv integer and is currently
% only supported for hypercubic graphs.
% ---------------------------------
% Format for NQS Modifier object with translation invariance:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nh + Alpha + 1. (computed here).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% Properties added with translation invariance:
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dWv.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Alpha = NQSObj.Alpha; % Hidden unit density, needs to be integer.

Cfg_vec = HilbertObj.FullCfgMod(Cfg); % Build the spin configuration vector.

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

if NQSObj.OptInds(1) == 1
    dLogp(1) = sum(Cfg_vec(:)); % Insert d/da.
end

% Accounting for shift structure of W matrix requires either construction
% of shifted Theta matrix or shifted Cfg vector - the latter is done here
for a=1:Alpha % Derivatives need to be computed by Alpha sector
    if NQSObj.OptInds(1+a) == 1
        dLogp(1+a) = sum(tanh(NQSObj.Theta((1:(2*Ntr))+(a-1)*2*Ntr))); % Insert d/db.
    end
    for v = 1:Nv
        PInd = 1 + Alpha + v + (a-1)*Nv;
        % For each layer labelled by a, find the indices of the associated translates.
        if NQSObj.OptInds(PInd) == 1
            for b = 1:2*Ntr
                TInd = b + (a-1)*2*Ntr; VInd = BondMap{1+mod(b-1,Ntr)}(1+mod(v-1,Ng)) + ...
                    Ng*mod(ceil(v/Ng)+ceil(b/Ntr),2);
                dLogp(PInd) = dLogp(PInd) + (Cfg_vec(VInd)*tanh(NQSObj.Theta(TInd))) + ...
                    (Cfg_vec(VInd)*tanh(NQSObj.Theta(TInd))); % Insert d/dW.
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
