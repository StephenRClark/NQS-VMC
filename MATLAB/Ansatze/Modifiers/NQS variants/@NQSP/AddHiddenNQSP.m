% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSP(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, W and Theta.
% ---------------------------------
% Format for NQSP Modifier:
% - NQSP.Nv = number of "visible" spins.
% - NQSP.Nh = number of "hidden" spins.
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

% Params requires field AlphaP, b, W, nphs, nmag.
AP = round(Params.AlphaP); % Require integer AlphaP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
A0 = NQSObj.Alpha; % Starting hidden unit density.
VOrder = NQSObj.VOrder; HOrder = NQSObj.HOrder;
OptInds = NQSObj.OptInds; % Optimisation indices, arranged [Re(p), Im(p)]

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Nsl = max(SLInds); Ntr = numel(BondMap); Ng = GraphObj.N;

bv = NQSObj.bv; Wm = NQSObj.Wm;

OptInds_a = OptInds(1:(Nsl*VOrder),:);
OptInds_b0 = OptInds(Nsl*VOrder+(1:(A0*HOrder)),:);
OptInds_W0 = OptInds(Nsl*VOrder+A0*HOrder+(1:(A0*Nv*HOrder*VOrder)),:);

if AP < 0
    if abs(AP) >= A0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AF = A0 + AP;
        bvF = bv(1:AF,:); OptInds_bF = OptInds_b0(1:(AF*HOrder),:);
        WmF = Wm(1:AF,:,:,:); OptInds_WF = OptInds_W0(1:(AF*Nv*HOrder*VOrder),:);
    end
else
    AF = A0 + AP;
    bvT = zeros(AP,HOrder); WmT = zeros(AP,Nv,HOrder,VOrder);
    for p = 1:AP
        for ho = 1:HOrder
            bvT(p,ho) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            for v = 1:Nv
                for vo = 1:VOrder
                    WmT(p,v,ho,vo) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
                end
            end
        end
    end
    OptInds_WT_R = real(permute(WmT,[4 3 2 1]))~=0; 
    OptInds_WT_I = imag(permute(WmT,[4 3 2 1]))~=0;
    OptInds_bT_R = real(bvT.')~=0; 
    OptInds_bT_I = imag(bvT.')~=0;
    OptInds_WT = [OptInds_WT_R(:), OptInds_WT_I(:)];
    OptInds_bT = [OptInds_bT_R(:), OptInds_bT_I(:)];
    bvF = cat(1,bv,bvT); WmF = cat(1,Wm,WmT);
    OptInds_bF = [OptInds_b0; OptInds_bT];
    OptInds_WF = [OptInds_W0; OptInds_WT];
end

Alpha = AF; NhF = Alpha*Ntr;

% Sort out optimisation indices
OptInds = [OptInds_a; OptInds_bF; OptInds_WF];

% Reassign all fields affected by Nh change.
NQSObj.Np = Nsl*VOrder + AF*HOrder + AF*Nv*HOrder*VOrder; NQSObj.Nh = AF; NQSObj.Theta = zeros(AF,1);
NQSObj.b = bvF; NQSObj.W = WmF; 
NQSObj.Alpha = Alpha; NQSObj.OptInds = OptInds;
end