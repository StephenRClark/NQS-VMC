% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSU(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will n-ecessitate changes in Nh, Np, b, W and Theta.
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

% Params requires field AlphaP, b, W, nphs, nmag.
AP = round(Params.AlphaP); % Require integer AlphaP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; Nmax = NQSObj.VDim-1; % Number of "visible" units and visible dimension.
A0 = NQSObj.Alpha; % Number of unique coupling sets.
OptInds = NQSObj.OptInds; % Optimisation indices, arranged [Re(p), Im(p)]

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

bv = NQSObj.bv; Wm = NQSObj.Wm; 

OptInds_a = OptInds(1:(Nsl*Nmax),:);
OptInds_b0 = OptInds(Nsl*Nmax+(1:A0),:);
OptInds_W0 = OptInds(Nsl*Nmax+A0+(1:(A0*Nv*Nmax)),:);

if AP < 0
    if abs(AP) >= A0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AF = A0 + AP;
        bvF = bv(1:AF); OptInds_bF = OptInds_b0(1:AF,:);     
        WmF = Wm(1:AF,:); OptInds_WF = OptInds_W0((1:(AF*Nv*Nmax)),:);
    end
else
    AF = A0 + AP;
    bvT = zeros(AP,1); WmT = zeros(AP,Nmax*Nv);
    for p = 1:AP
        bvT(p) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nmax*Nv
            WmT(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    OptInds_WT_R = (real(WmT.')~=0); OptInds_WT_I = (imag(WmT.')~=0);
    bvF = [bv; bvT]; OptInds_bT = [(real(bvT)~=0), (imag(bvT)~=0)];
    WmF = [Wm; WmT]; OptInds_WT = [OptInds_WT_R(:), OptInds_WT_I(:)];
    OptInds_bF = [OptInds_b0; OptInds_bT];
    OptInds_WF = [OptInds_W0; OptInds_WT];
end

Alpha = AF; NhF = Alpha*Ntr;

% Sort out optimisation indices
OptInds = [OptInds_a; OptInds_bF; OptInds_WF];

% Reassign all fields affected by Nh change.
NQSObj.Np = Nmax*Nsl + AF + AF*Nmax*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.bv = bvF; NQSObj.Wm = WmF; 
NQSObj.Alpha = Alpha; NQSObj.OptInds = OptInds;

% Repackage the parameters in the necessary form.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for bd = 1:Ntr
        for n = 1:Nv
            if BondMap{bd}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{bd}(1+mod(n-1,Ng));
                for v = 1:Nmax
                    NInd = v + Nmax*(n-1); WInd = v + Nmax*(VInd-1);
                    % Account for enlarged lattices where Nv = Ns x Ng.
                    NQSObj.W(bd+(al-1)*Ntr,WInd) = NQSObj.Wm(al,NInd);
                end
            end
        end
    end
end

end