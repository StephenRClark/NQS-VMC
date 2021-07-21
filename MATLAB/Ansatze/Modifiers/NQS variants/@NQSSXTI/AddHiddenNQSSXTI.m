% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSSXTI(NQSObj,Params)
% This function adds round(NhP/Nv) new hidden units to an existing NQSObj
% (removes if negative). This will necessitate changes in Nh, Np, b, bti,
% W, Wv, Alpha and Theta.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% Properties added with translation invariance:
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------

% Params requires field NhP or AlphaP, b, W, nphs, nmag.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Alpha0 = NQSObj.Alpha; % Starting hidden unit density.
bti = NQSObj.bti; Bti = NQSObj.Bti; Xv = NQSObj.Xv; Wv = NQSObj.Wv;
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

if isfield(Params,'AlphaP')
    AlphaP = Params.AlphaP;
elseif isfield(Params,'NhP')
    AlphaP = round(NhP/Nv);
end

if isfield(Params,'B') == 0
    Params.B = Params.b;
end
if isfield(Params,'X') == 0
    Params.X = Params.W;
end

if AlphaP < 0
    if abs(AlphaP) >= Alpha0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AlphaF = Alpha0 + AlphaP;
        btiF = bti(1:AlphaF); XvF = Xv(1:AlphaF,:); WvF = Wv(1:AlphaF,:);
    end
else
    AlphaF = Alpha0 + AlphaP;
    btiF = zeros(AlphaP,1); XvF = zeros(AlphaP,Nv); WvF = zeros(AlphaP,Nv);
    for a = 1:AlphaP
        btiF(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BtiF(a) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WvF(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            XvF(a,v) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);            
        end
    end
    btiF = [bti; btiF]; BtiF = [Bti; BtiF]; XvF = [Xv; XvF]; WvF = [Wv; WvF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = 2 + 2*AlphaF + 2*AlphaF*Nv; NQSObj.Nh = AlphaF*Ntr; NQSObj.Theta = zeros(AlphaF*Ntr,1);
NQSObj.bti = btiF; NQSObj.Bti = BtiF; NQSObj.Xv = XvF; NQSObj.Wv = WvF; 
NQSObj.OptInds = ones(NQSObj.Np,1); NQSObj.Alpha = AlphaF;

% Constructing shift invariant W matrix.
for a = 1:AlphaF
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    NQSObj.B((1:Ntr)+(a-1)*Ntr) = NQSObj.Bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.X(b+(a-1)*Ntr,VInd) = NQSObj.Xv(a,n);
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wv(a,n);
            end
        end
    end
end

end