% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSS1(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, B, W and Theta.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Params requires field AlphaP, b, W, nphs, nmag.
AP = round(Params.AlphaP); % Require integer AlphaP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
A0 = NQSObj.Alpha; % Density of "hidden" units.
OptInds = NQSObj.OptInds; % Optimisation indices, arranged [Re(p), Im(p)]

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

OptInds_aA = OptInds(1:(2*Nsl),:);
OptInds_b0 = OptInds(2*Nsl+(1:A0),:);
OptInds_w0 = OptInds(2*Nsl+A0+(1:(A0*Nv)),:);
OptInds_W0 = OptInds(2*Nsl+A0*(1+Nv)+(1:(A0*Nv)),:);

if isfield(Params,'w') == 0
    Params.w = Params.W;
end

bv = NQSObj.bv; Wm = NQSObj.Wm; wm = NQSObj.wm;

if AP < 0 % Remove hidden units.
    if abs(AP) >= A0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AF = A0 + AP;
        bvF = bv(1:AF); OptInds_bvF = OptInds_b0(1:AF,:);  
        wmF = wm(1:AF,:); OptInds_wmF = OptInds_w0((1:(AF*Nv)),:);
        WmF = Wm(1:AF,:); OptInds_WmF = OptInds_W0((1:(AF*Nv)),:);
    end
else
    AF = A0 + AP;
    bvT = zeros(AP,1); WmT = zeros(AP,Nv); wmT = zeros(AP,Nv);
    for a = 1:AP
        bvT(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            wmT(a,v) = Params.w * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            WmT(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    OptInds_wmT_R = (real(wmT.')~=0); OptInds_wmT_I = (imag(wmT.')~=0);
    wmF = [wm; wmT]; OptInds_wmT = [OptInds_wmT_R(:), OptInds_wmT_I(:)];
    OptInds_WmT_R = (real(WmT.')~=0); OptInds_WmT_I = (imag(WmT.')~=0);
    WmF = [Wm; WmT]; OptInds_WT = [OptInds_WmT_R(:), OptInds_WmT_I(:)];
    bvF = [bv; bvT]; OptInds_bT = [(real(bvT)~=0), (imag(bvT)~=0)];
    OptInds_bvF = [OptInds_b0; OptInds_bT];
    OptInds_wmF = [OptInds_w0; OptInds_wmT];
    OptInds_WmF = [OptInds_W0; OptInds_WT];
end

Alpha = AF; NhF = Alpha*Ntr;

% Sort out optimisation indices
OptInds = [OptInds_aA; OptInds_bvF; OptInds_wmF; OptInds_WmF];

% Reassign all fields affected by Nh change.
NQSObj.Np = 2*Nsl + AF + 2*AF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.bv = bvF; NQSObj.wm = wmF; NQSObj.Wm = WmF; NQSObj.OptInds = OptInds;

% Repackage the ati, bti and Wv to usual NQS form.

% Constructing shift invariant W matrix.
for a = 1:A0
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bv(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.                
                NQSObj.w(b+(a-1)*Ntr,VInd) = NQSObj.wm(a,n);
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wm(a,n);
            end
        end
    end
end

end