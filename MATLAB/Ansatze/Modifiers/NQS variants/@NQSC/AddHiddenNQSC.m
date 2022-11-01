% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSC(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, B, W and Theta.
% ---------------------------------
% Format for NQSC Modifier:
% - NQSC.Nv = number of "visible" units.
% - NQSC.Nh = number of "hidden" units.
% - NQSC.Np = number of parameters in the ansatz = 2*Alpha + 2*Alpha*Nv + 2*Nsl.
% - NQSC.a = (Nv x 1) vector - visible site bias.
% - NQSC.av = (Nsl x 1) vector - visible bias parameters.
% - NQSC.A = (Nv x 1) vector - visible site square bias.
% - NQSC.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSC.b = (Nh x 1) vector - hidden site bias.
% - NQSC.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSC.B = (Nh x 1) vector- hidden site square bias.
% - NQSC.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSC.w = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSC.wm = (Alpha x Nv) matrix - coupling parameters
% - NQSC.W = (Nh x Nv) matrix - hidden-square-visible coupling terms.
% - NQSC.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSC.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSC.HDim = dimension of the hidden units.
% - NQSC.Theta = (Nh x 1) vector - effective angles.
% - NQSC.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
    
% Params requires field AlphaP, b, W, nphs, nmag.
AP = round(Params.AlphaP); % Require integer AlphaP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
A0 = NQSObj.Alpha; % Starting hidden unit density.
OptInds = NQSObj.OptInds; % Optimisation indices, arranged [Re(p), Im(p)]

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

bv = NQSObj.bv; Bv = NQSObj.Bv; wm = NQSObj.wm; Wm = NQSObj.Wm;

OptInds_aA = OptInds(1:(2*Nsl),:);
OptInds_b0 = OptInds(2*Nsl+(1:A0),:);
OptInds_B0 = OptInds(2*Nsl+A0+(1:A0),:);
OptInds_w0 = OptInds(2*Nsl+2*A0+(1:(A0*Nv)),:);
OptInds_W0 = OptInds(2*Nsl+2*A0+A0*Nv+(1:(A0*Nv)),:);

if ~isfield(Params,'B')
    Params.B = Params.b;
end
if ~isfield(Params,'w')
    Params.w = Params.W;
end

if AP < 0 % Remove hidden units.
    if abs(AP) >= A0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AF = A0 + AP;
        bvF = bv(1:AF); OptInds_bF = OptInds_b0(1:AF,:);  
        BvF = Bv(1:AF); OptInds_BF = OptInds_B0(1:AF,:);        
        wmF = wm(1:AF,:); OptInds_wF = OptInds_w0((1:(AF*Nv)),:);
        WmF = Wm(1:AF,:); OptInds_WF = OptInds_W0((1:(AF*Nv)),:);
    end
else
    AF = A0 + AP;
    bvT = zeros(AP,1); BvT = zeros(AP,1); wmT = zeros(AP,Nv); WmT = zeros(AP,Nv);
    for a = 1:AP
        bvT(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BvT(a) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            wmT(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            WmT(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    OptInds_wT_R = (real(wmT.')~=0); OptInds_wT_I = (imag(wmT.')~=0);
    OptInds_WT_R = (real(WmT.')~=0); OptInds_WT_I = (imag(WmT.')~=0);
    bvF = [bv; bvT]; OptInds_bT = [(real(bvT)~=0), (imag(bvT)~=0)];
    BvF = [Bv; BvT]; OptInds_BT = [(real(BvT)~=0), (imag(BvT)~=0)];
    wmF = [wm; wmT]; OptInds_wT = [OptInds_wT_R(:), OptInds_wT_I(:)];
    WmF = [Wm; WmT]; OptInds_WT = [OptInds_WT_R(:), OptInds_WT_I(:)];
    OptInds_bF = [OptInds_b0; OptInds_bT];
    OptInds_BF = [OptInds_B0; OptInds_BT];
    OptInds_wF = [OptInds_w0; OptInds_wT];
    OptInds_WF = [OptInds_W0; OptInds_WT];
end

Alpha = AF; NhF = Alpha*Ntr;

% Sort out optimisation indices
OptInds = [OptInds_aA; OptInds_bF; OptInds_BF; OptInds_wF; OptInds_WF];

% Reassign all fields affected by Nh change.
NQSObj.Np = 2*Nsl + 2*AF + 2*AF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.bv = bvF; NQSObj.Bv = BvF; NQSObj.wm = wmF; NQSObj.Wm = WmF; 
NQSObj.Alpha = Alpha; NQSObj.OptInds = OptInds;

% Constructing shift invariant W matrix.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    NQSObj.B((1:Ntr)+(al-1)*Ntr) = NQSObj.Bv(al);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.w(b+(al-1)*Ntr,VInd) = NQSObj.wm(al,n);
                NQSObj.W(b+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
            end
        end
    end
end

end