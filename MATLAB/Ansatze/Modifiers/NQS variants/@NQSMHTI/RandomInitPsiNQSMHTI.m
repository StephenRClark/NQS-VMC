% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSMHTI(NQSObj,Params)
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
% - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
% - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
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
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% NB: Translation invariance requires the desired translation vectors to be
% specified in the Graph provided.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Alpha = round(Nh/Ng); NQSObj.Alpha = Alpha; % Hidden unit density, needs to be integer.
NQSObj.Np = 2*Nv*Alpha + 2*Alpha + 2; % The number of variational parameters.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.

Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

Nh = Ntr * Alpha; NQSObj.Nh = Nh; % Determine and reassign Nh required.

% Initialise the storage:

NQSObj.BH = zeros(Nh,1);
NQSObj.BHti = zeros(Alpha,1);
NQSObj.BM = zeros(Nh,1);
NQSObj.BMti = zeros(Alpha,1);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Wv = zeros(Alpha,Nv);
NQSObj.X = zeros(Nh,Nv);
NQSObj.Xv = zeros(Alpha,Nv);
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

NQSObj.ati = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
NQSObj.Ati = (Params.A + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.A~=0);

for a = 1:Alpha
    NQSObj.BHti(a) = (Params.BH + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.BH~=0);
    NQSObj.BMti(a) = (Params.BM + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.BM~=0);
    for n = 1:Nv
        NQSObj.Wv(a,n) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
        NQSObj.Xv(a,n) = (Params.X + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.X~=0);
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
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wv(a,n);
                NQSObj.X(b+(a-1)*Ntr,VInd) = NQSObj.Xv(a,n);
            end
        end
    end
end

NQSObj.OptInds = ones(NQSObj.Np,1);
end