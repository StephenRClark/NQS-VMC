% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSMHAI(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% Ansatz is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% NB: Translation invariance here assumes Nh/Nv integer.
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
% NB: Translation invariance requires the desired translation vectors to be
% specified in the Graph provided.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Alpha = round(Nh/Ng); NQSObj.Alpha = Alpha; % Hidden unit density, needs to be integer.
NQSObj.Np = 4*Nv*Alpha + 2*Alpha + 2; % The number of variational parameters.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.

Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

Nh = Ntr * Alpha; NQSObj.Nh = Nh; % Determine and reassign Nh required.

% Initialise the storage:

NQSObj.BH = zeros(Nh,1);
NQSObj.BHti = zeros(Alpha,1);
NQSObj.BM = zeros(Nh,1);
NQSObj.BMti = zeros(Alpha,1);
NQSObj.WH = zeros(Nh,Nv);
NQSObj.WM = zeros(Nh,Nv);
NQSObj.WHv = zeros(Alpha,Nv);
NQSObj.WMv = zeros(Alpha,Nv);
NQSObj.XH = zeros(Nh,Nv);
NQSObj.XM = zeros(Nh,Nv);
NQSObj.XHv = zeros(Alpha,Nv);
NQSObj.XMv = zeros(Alpha,Nv);
NQSObj.ThetaH = zeros(Nh,1);
NQSObj.ThetaM = zeros(Nh,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
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

NQSObj.ati = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
NQSObj.Ati = Params.A * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);

for a = 1:Alpha
    NQSObj.BHti(a) = Params.BH * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    NQSObj.BMti(a) = Params.BM * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    for n = 1:Nv
        NQSObj.WHv(a,n) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        NQSObj.WMv(a,n) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        NQSObj.XHv(a,n) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        NQSObj.XMv(a,n) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    end
end

% Repackage the ati, bti and Wv to usual NQS form.
NQSObj.a = NQSObj.ati * ones(Nv,1);
NQSObj.A = NQSObj.Ati * ones(Nv,1);

% Constructing shift invariant W matrix.
for a = 1:Alpha
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

NQSObj.OptInds = ones(NQSObj.Np,1);