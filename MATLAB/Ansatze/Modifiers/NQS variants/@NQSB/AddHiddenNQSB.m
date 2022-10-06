% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSB(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, B, W and Theta.
% ---------------------------------
% Format for NQSB Modifier:
% - NQSB.Nv = number of "visible" units.
% - NQSB.Nh = number of "hidden" units.
% - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSB.a = (Nv x 1) vector - visible site bias.
% - NQSB.av = (Nsl x 1) vector - visible bias parameters.
% - NQSB.A = (Nv x 1) vector - visible site square bias.
% - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSB.b = (Nh x 1) vector - hidden site bias.
% - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSB.B = (Nh x 1) vector- hidden site square bias.
% - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSB.HDim = dimension of the hidden units.
% - NQSB.Theta = (Nh x 1) vector - effective angles.
% - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
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

bv = NQSObj.bv; Bv = NQSObj.Bv; Wm = NQSObj.Wm;

OptInds_a = OptInds(1:(2*Nsl),:);
OptInds_b0 = OptInds(2*Nsl+(1:A0),:);
OptInds_B0 = OptInds(2*Nsl+A0+(1:A0),:);
OptInds_W0 = OptInds(2*Nsl+2*A0+(1:(A0*Nv)),:);

if isfield(Params,'B') == 0
    Params.B = Params.b;
end

if AP < 0 % Remove hidden units.
    if abs(AP) >= A0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AF = A0 + AP;
        bvF = bv(1:AF); OptInds_bF = OptInds_b0(1:AF,:);  
        BvF = Bv(1:AF); OptInds_BF = OptInds_B0(1:AF,:);        
        WmF = Wm(1:AF,:); OptInds_WF = OptInds_W0((1:(AF*Nv)),:);
    end
else
    AF = A0 + AP;
    bvT = zeros(AP,1); BvT = zeros(AP,1); WmT = zeros(AP,Nv);
    for a = 1:AP
        bvT(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BvT(a) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WmT(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    OptInds_WT_R = (real(WmT.')~=0); OptInds_WT_I = (imag(WmT.')~=0);
    bvF = [bv; bvT]; OptInds_bT = [(real(bvT)~=0), (imag(bvT)~=0)];
    BvF = [Bv; BvT]; OptInds_BT = [(real(BvT)~=0), (imag(BvT)~=0)];
    WmF = [Wm; WmT]; OptInds_WT = [OptInds_WT_R(:), OptInds_WT_I(:)];
    OptInds_bF = [OptInds_b0; OptInds_bT];
    OptInds_BF = [OptInds_B0; OptInds_BT];
    OptInds_WF = [OptInds_W0; OptInds_WT];
end

Alpha = AF; NhF = Alpha*Ntr;

% Sort out optimisation indices
OptInds = [OptInds_a; OptInds_bF; OptInds_BF; OptInds_WF];

% Reassign all fields affected by Nh change.
NQSObj.Np = 2*Nsl + 2*AF + AF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.bv = bvF; NQSObj.Bv = BvF; NQSObj.Wm = WmF; 
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
                NQSObj.W(b+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
            end
        end
    end
end

end