% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSU(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" units.
% - NQSU.Nh = number of "hidden" units.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nmax*Nsl x 1) for d/da.
% Arranged [sl, vd], [sl, vd+1], ... , [sl+1, vd], ...
% - (Alpha x 1) for d/db.
% - (Alpha*Nv*Nmax x 1) for d/dW.
% Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; Nmax = NQSObj.VDim-1; % Number of "visible" units and visible dimension.
% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;

Alpha = NQSObj.Alpha; % Number of unique coupling sets.
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
UMat = zeros(Nmax,Nv);

dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.

for v = 1:Nmax
    UMat(v,:) = (Cfg_vec.' == NQSObj.VList(v+1));
    for s = 1:Nsl
        PInd = v + (s-1)*Nmax;
        if sum(NQSObj.OptInds(PInd,:))~=0
            % Use sublattice indices to determine which sites contribute
            dLogp(PInd) = sum(UMat(PInd,GraphObj.SLInds==s));
        end
    end
end

dTheta = tanh(NQSObj.Theta);
for a = 1:Alpha
    bInd = Nmax*Nsl + a;
    if sum(NQSObj.OptInds(bInd,:)) ~= 0
        dLogp(bInd) = sum(dTheta((1:Ntr)+(a-1)*Ntr));
    end
    for n = 1:Nv
        for v = 1:Nmax
            PInd = Nmax*Nsl + Alpha + Nmax*((n-1)+(a-1)*Nv) + v;
            if sum(NQSObj.OptInds(PInd,:)) ~= 0
                for b = 1:numel(BondMap)
                    HInd = b + (a-1)*Ntr; VInd = BondMap{b}(n);
                    dLogp(PInd) = dLogp(PInd) + UMat(v,VInd)*dTheta(HInd);
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
end