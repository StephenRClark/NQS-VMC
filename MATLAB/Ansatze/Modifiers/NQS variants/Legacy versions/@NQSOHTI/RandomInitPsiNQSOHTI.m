% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSOHTI(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
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

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = NQSObj.Nh; % Number of "hidden" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Alpha = round(Nh/Ng); NQSObj.Alpha = Alpha; % Hidden unit density, needs to be integer.
NQSObj.Np = Nv*Alpha + Alpha + Nv/Ng; % The number of variational parameters.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Nh = Ntr * Alpha; NQSObj.Nh = Nh; % Determine and reassign Nh required.

NQSObj.Np = VDim + Alpha + (VDim*Nv * Alpha); % The number of variational parameters.

% Initialise the storage:
NQSObj.a = zeros(VDim*Nv,1);
NQSObj.b = zeros(Nh,1);
NQSObj.W = zeros(Nh,VDim*Nv);
NQSObj.Theta = zeros(Nh,1);
NQSObj.ati = zeros(VDim,1);
NQSObj.bti = zeros(Alpha,1);
NQSObj.Wv = zeros(Alpha,VDim*Nv);

for v = 1:VDim
    NQSObj.ati(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
end
for a = 1:Alpha
    NQSObj.bti(a) = (Params.b + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
end
for a = 1:Alpha
    for v = 1:VDim*Nv
        NQSObj.Wv(a,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
    end
end

NQSObj.OptInds = ones(NQSObj.Np,1);

% Repackage the parameters in the necessary form.
NQSObj.a = reshape(NQSObj.ati * ones(1,Nv),VDim*Nv,1);
for a = 1:Alpha
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