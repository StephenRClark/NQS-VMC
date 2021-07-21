% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSNHTI(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% Ansatz is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% NB: Translation invariance here assumes Nh/Nv integer.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nh + 2*Alpha + 2.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.B = (Nh x 1) vector - hidden site square bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% Properties added with translation invariance:
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
% ---------------------------------
% NB: Translation invariance requires the desired translation vectors to be
% specified in the Graph provided.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Alpha = round(Nh/Ng); NQSObj.Alpha = Alpha; % Hidden unit density, needs to be integer.
NQSObj.Np = Nv*Alpha + 2*Alpha + 2; % The number of variational parameters.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.

Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

Nh = Ntr * Alpha; NQSObj.Nh = Nh; % Determine and reassign Nh required.

% Initialise the storage:

NQSObj.b = zeros(Nh,1);
NQSObj.bti = zeros(Alpha,1);
NQSObj.B = zeros(Nh,1);
NQSObj.Bti = zeros(Alpha,1);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Wv = zeros(Alpha,Nv);
NQSObj.Theta = zeros(Nh,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
end
if isfield(Params,'B') == 0
    Params.B = Params.b;
end

NQSObj.ati = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
NQSObj.Ati = (Params.A + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.A~=0);
for a = 1:Alpha
    NQSObj.bti(a) = (Params.b + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
    NQSObj.Bti(a) = (Params.B + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.B~=0);
    for n = 1:Nv
        NQSObj.Wv(a,n) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
    end
end

% Repackage the ati, bti and Wv to usual NQS form.
NQSObj.a = NQSObj.ati * ones(Nv,1);
NQSObj.A = NQSObj.Ati * ones(Nv,1);

% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    NQSObj.B((1:Ntr)+(a-1)*Ntr) = NQSObj.Bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wv(a,n);
            end
        end
    end
end

NQSObj.OptInds = ones(NQSObj.Np,1);
end