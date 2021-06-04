% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSOHTI(NQSObj,Params)
% This function adds round(NhP/Nv) new hidden units to an existing NQSObj
% (removes if negative). This will necessitate changes in Nh, Np, b, bti,
% W, Wv, Alpha and Theta.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim + Alpha + (VDim*Nv * Alpha).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% Properties added with translation invariance:
% - NQS.ati = (VDim x 1) vector - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x VDim*Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------

% Params requires field NhP or AlphaP, b, c, nphs, nmag.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
Alpha0 = NQSObj.Alpha; % Starting hidden unit density.
bti = NQSObj.bti; Wv = NQSObj.Wv;
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

if isfield(Params,'AlphaP')
    AlphaP = Params.AlphaP;
elseif isfield(Params,'NhP')
    AlphaP = round(NhP/Nv);
end

if AlphaP < 0
    if abs(AlphaP) >= Alpha0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AlphaF = Alpha0 + AlphaP;
        btiF = bti(1:AlphaF); WvF = Wv(1:AlphaF,:);
    end
else
    AlphaF = Alpha0 + AlphaP;
    btiF = zeros(AlphaP,1); WvF = zeros(AlphaP,VDim*Nv);
    for a = 1:AlphaP
        btiF(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:VDim*Nv
            WvF(a,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    btiF = [bti; btiF]; WvF = [Wv; WvF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = VDim + AlphaF + AlphaF*VDim*Nv; NQSObj.Nh = AlphaF*Ntr; NQSObj.Theta = zeros(AlphaF*Ntr,1);
NQSObj.bti = btiF; NQSObj.Wv = WvF; NQSObj.OptInds = ones(NQSObj.Np,1); NQSObj.Alpha = AlphaF;

% Constructing shift invariant W matrix.
for a = 1:AlphaF
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng));
                for v = 1:VDim
                    NInd = v + VDim*(n-1); WInd = v + VDim*(VInd-1);
                    % Account for enlarged lattices where Nv = Ns x Ng.
                    NQSObj.W(b+(a-1)*Ntr,WInd) = NQSObj.Wv(a,NInd);
                end
            end
        end
    end
end

end