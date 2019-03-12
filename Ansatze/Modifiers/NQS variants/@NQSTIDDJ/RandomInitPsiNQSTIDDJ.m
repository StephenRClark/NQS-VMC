% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSTIDDJ(NQSObj,GraphObj,Params)
% This function populates random initial NQS ansatz structure. The input
% Ansatz is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% NB: Translation invariance here assumes Nh/Nv integer, and spin symmetry
% is imposed. Spin symmetry results in effectively double the hidden unit
% number.
% ---------------------------------
% Format for NQS Modifier object with translation invariance:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv*Alpha + Alpha + Nv/Ng. (computed here).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% Properties added with translation invariance:
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------
% NB: Translation invariance requires the desired translation vectors to be
% specified in the Graph provided.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.
Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Alpha = round(Nh/Ng); % Hidden unit density, needs to be integer.
NQSObj.Alpha = Alpha; % Reassign Alpha, as Nh will be changed.
NQSObj.Np = (Ng+1)*Alpha + Alpha + 1; % The number of variational parameters.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.

Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

Nh = Ntr * Alpha; NQSObj.Nh = 2*Nh; % Determine and reassign Nh required.

% Initialise the storage:
NQSObj.b = zeros(2*Nh,1);
NQSObj.bti = zeros(Alpha,1);
NQSObj.W = zeros(2*Nh,Nv);
NQSObj.Wv = zeros(Alpha,Ng+1);
NQSObj.Theta = zeros(2*Nh,1);

NQSObj.ati = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);

for a = 1:Alpha
    NQSObj.bti(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    for n = 1:(Ng+1)
        NQSObj.Wv(a,n) = Params.c * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    end
end

% Repackage the ati, bti and Wv to usual NQS form.
NQSObj.a = NQSObj.ati * ones(Nv,1);

% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Ng
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng));
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wv(a,n);
                if n == 1
                    NQSObj.W(b+(a-1)*Ntr,VInd+Ng) = NQSObj.Wv(a,Ng+1);
                else
                    NQSObj.W(b+(a-1)*Ntr,VInd+Ng)= NQSObj.Wv(a,n);
                end
            end
        end
    end
end

NQSObj.b((1:Nh)+Nh) = NQSObj.b(1:Nh);

NQSObj.W((1:Nh)+Nh,:) = [NQSObj.W(1:Nh,(1:Ng)+Ng) NQSObj.W(1:Nh,1:Ng)];

NQSObj.OptInds = ones(NQSObj.Np,1);