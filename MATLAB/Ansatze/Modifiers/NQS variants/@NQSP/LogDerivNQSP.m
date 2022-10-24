% --- General NQS logarithmic derivative function ---

function dLogp = LogDerivNQSP(NQSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for NQSP Modifier:
% - NQSP.Nv = number of "visible" units.
% - NQSP.Nh = number of "hidden" units.
% - NQSP.Np = number of parameters in the ansatz = (Nsl x VOrder) + (Alpha x
% HOrder) + (Nv x VOrder)(Alpha x HOrder)
% - NQSP.VDim = dimension of the visible units.
% - NQSP.HDim = dimension of the hidden units.
% - NQSP.VOrder = highest power of visible unit interactions. Max value VDim-1.
% - NQSP.HOrder = highest power of hidden unit interactions. Max value HDim-1.
% - NQSP.a = (Nv x VOrder) matrix - visible site biases.
% - NQSP.av = (Nsl x VOrder) matrix - visible bias parameters
% - NQSP.b = (Nh x HOrder) matrix - hidden site bias.
% - NQSP.bv = (Alpha x HOrder) matrix - hidden bias parameters.
% - NQSP.W = (Nh x Nv x HOrder x VOrder) array - hidden-visible coupling terms.
% - NQSP.Wm = (Alpha x Nv x HOrder x VOrder) array - hidden-visible coupling parameters
% - NQSP.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSP.Theta = (Nh x HOrder) matrix - effective angles by hidden order.
% - NQSP.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSP.Rescale = flag for visible unit rescaling to [0 1] interval.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x VOrder) x 1 for d/da. Group by Sublattice > Visible order
% > [sl,vo], [sl, vo+1] ... [sl+1, vo]
% - (Alpha x HOrder) x 1 for d/db. Group by Alpha > Hidden order
% > [al, ho], [al, ho+1] ... [al+1, ho]
% - (Alpha x Nv) x (HOrder x VOrder) for d/dW. Group by Alpha > Position > Hidden order > Visible order
% > [al,v,ho,vo], [al,v,ho,vo+1] ... [al,v,ho+1,vo] ... [al,v+1,ho,vo] ... [al+1,v,ho,vo]
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.
HOrder = NQSObj.HOrder; HDim = NQSObj.HDim;
VOrder = NQSObj.VOrder; VDim = NQSObj.VDim;
Rescale = NQSObj.Rescale;

GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;

Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

Cfg_vec = NQSObj.FullCfg(Cfg)*((VDim-1)^-Rescale); % Build the spin configuration vector.

NQSObj.VisVec = Cfg_vec; % Ensure VisVec is properly assigned.
dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
VisPow = Cfg_vec(:) .^ (1:VOrder); % Nv x VOrder
HidPow = ((0:(HDim-1)).'*((HDim-1)^(-Rescale))) .^ (1:HOrder); % HDim x HOrder
HidArray = reshape(HidPow.',1,HOrder,HDim);
dTheta = (sum(exp(NQSObj.Theta.*HidArray).*HidArray,3))./sum(exp(NQSObj.Theta.*HidArray),3); % Nh x HOrder
for s = 1:Nsl
    for vo = 1:VOrder
        PInd = vo + (s-1)*VOrder;
        if sum(NQSObj.OptInds(PInd,:))~=0
            dLogp(PInd) = sum(VisPow(SLInds==s,vo)); % Insert d/da.
        end
    end
end
for al = 1:Alpha
    for ho = 1:HOrder
        PInd = Nsl*VOrder + ho + (al-1)*HOrder;
        if sum(NQSObj.OptInds(PInd,:))~=0
            dLogp(PInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr,ho));
        end
        for v = 1:Nv
            for vo = 1:VOrder
                PInd = Nsl*VOrder + Alpha*HOrder + ...
                    vo + VOrder(ho-1 + HOrder*(v-1 + (al-1)*Nv));
                if sum(NQSObj.OptInds(PInd,:))~=0
                    for bd = 1:numel(BondMap)
                        HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                        if VInd ~= 0
                            dLogp(PInd) = dLogp(PInd) + VisPow(VInd,vo)*dTheta(HInd,ho); % Insert d/dW
                        end
                    end
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
end