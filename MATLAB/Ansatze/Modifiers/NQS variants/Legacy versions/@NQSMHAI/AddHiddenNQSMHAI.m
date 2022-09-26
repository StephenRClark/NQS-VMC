% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSMHAI(NQSObj,Params)
% This function adds round(NhP/Nv) new hidden units to an existing NQSObj
% (removes if negative). This will necessitate changes in Nh, Np, b, bti,
% W, Wv, Alpha and Theta.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.WH = (Nh x Nv) matrix - hidden-visible HH coupling terms.
% - NQS.WM = (Nh x Nv) matrix - hidden-visible MM coupling terms.
% - NQS.XH = (Nh x Nv) matrix - hidden-visible HM coupling terms.
% - NQS.XM = (Nh x Nv) matrix - hidden-visible MH coupling terms.
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
% - NQS.WHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.WMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------

% Params requires field NhP or AlphaP, b, W, nphs, nmag.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Alpha0 = NQSObj.Alpha; % Starting hidden unit density.
WHv = NQSObj.WHv; WMv = NQSObj.WMv; BHti = NQSObj.BHti;
XHv = NQSObj.XHv; XMv = NQSObj.XMv; BMti = NQSObj.BMti; 
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

if isfield(Params,'AlphaP')
    AlphaP = Params.AlphaP;
elseif isfield(Params,'NhP')
    AlphaP = round(NhP/Nv);
end

if (isfield(Params,'BH') == 0) && isfield(Params,'BM')
    Params.BH = Params.BM;
elseif (isfield(Params,'BM') == 0) && isfield(Params,'BH')
    Params.BM = Params.BH;
elseif (isfield(Params,'BH') == 0) && (isfield(Params,'BM') == 0) && isfield(Params,'b')
    Params.BH = Params.b; Params.BM = Params.b;
end
if isfield(Params,'X') == 0
    Params.X = Params.W;
end

if AlphaP < 0
    if abs(AlphaP) >= Alpha0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AlphaF = Alpha0 + AlphaP;        
        WHvF = WHv(1:AlphaF,:); WMvF = WMv(1:AlphaF,:); BHtiF = BHti(1:AlphaF); 
        XHvF = XHv(1:AlphaF,:); XMvF = XMv(1:AlphaF,:); BMtiF = BMti(1:AlphaF); 
    end
else
    AlphaF = Alpha0 + AlphaP;    
    WHvF = zeros(AlphaP,Nv); WMvF = zeros(AlphaP,Nv); BHtiF = zeros(AlphaP,1); 
    XHvF = zeros(AlphaP,Nv); XMvF = zeros(AlphaP,Nv); BMtiF = zeros(AlphaP,1); 
    for a = 1:AlphaP
        BHtiF(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BMtiF(a) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WHvF(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            WMvF(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            XHvF(a,v) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            XMvF(a,v) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end    
    WHvF = [WHv; WHvF]; WMvF = [WMv; WMvF]; BHtiF = [BHti; BHtiF];
    XHvF = [XHv; XHvF]; XMvF = [XMv; XMvF]; BMtiF = [BMti; BMtiF]; 
end

% Reassign all fields affected by Nh change.
NQSObj.Np = 2 + 2*AlphaF + 4*AlphaF*Nv; NQSObj.Nh = AlphaF*Ntr;
NQSObj.ThetaH = zeros(AlphaF*Ntr,1); NQSObj.ThetaM = zeros(AlphaF*Ntr,1);
NQSObj.WHv = WHvF; NQSObj.WMv = WMvF; NQSObj.BHti = BHtiF;
NQSObj.XHv = XHvF; NQSObj.XMv = XMvF; NQSObj.BMti = BMtiF;
NQSObj.OptInds = ones(NQSObj.Np,1); NQSObj.Alpha = AlphaF;

% Constructing shift invariant W matrix.
for a = 1:AlphaF
    NQSObj.BH((1:Ntr)+(a-1)*Ntr) = NQSObj.BHti(a);
    NQSObj.BM((1:Ntr)+(a-1)*Ntr) = NQSObj.BMti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.WH(b+(a-1)*Ntr,VInd) = NQSObj.WHv(a,n);
                NQSObj.WM(b+(a-1)*Ntr,VInd) = NQSObj.WMv(a,n);
                NQSObj.XH(b+(a-1)*Ntr,VInd) = NQSObj.XHv(a,n);
                NQSObj.XM(b+(a-1)*Ntr,VInd) = NQSObj.XMv(a,n);
            end
        end
    end
end