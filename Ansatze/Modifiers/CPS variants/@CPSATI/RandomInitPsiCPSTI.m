% --- General CPS wave function random initialisation function ---

% - Original by X Fang, updated by M Pei.

function [CPSObj] = RandomInitPsiCPSTI(CPSObj,Params)
% This function populates random initial CPS ansatz structure. The input
% Ansatz is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for CPS Modifier object:
% - CPS.Nv = number of "visible" spins.
% - CPS.Nh = number of "hidden" spins.
% - CPS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPS.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPS.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPS.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPS.HDim = 3 - this version features fixed hidden unit dimension.
% - CPS.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPS.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPS.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPS.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPS.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% Properties added with translation invariance:
% - CPS.Alpha = number of distinct hidden unit sets.
% - CPS.ati = (1 x (VDim-1)) vector - reduced parameter set for TI.
% - CPS.bti = (Alpha x (HDim-1)) matrix - reduced parameter set for TI.
% - CPS.Wti = ((VDim-1) x (HDim-1) x Nv x Alpha) array - reduced parameter set for TI.
% ---------------------------------

Nv = CPSObj.Nv; VDim = CPSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = CPSObj.Nh; HDim = CPSObj.HDim; % Number of "hidden" spins and hidden dimension.
GraphObj = CPSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
Alpha = round(Nh/Ng); CPSObj.Alpha = Alpha; % Hidden unit density, needs to be integer.
CPSObj.Np = Nv*Alpha + Alpha + Nv/Ng; % The number of variational parameters.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Nh = Ntr * Alpha; CPSObj.Nh = Nh; % Determine and reassign Nh required.

CPSObj.Np = (VDim-1) + (HDim-1) * Alpha + ((HDim-1)*(VDim-1) * (Nv * Alpha)); % The number of variational parameters.
CPSObj.a = zeros(Nv,VDim-1);
CPSObj.b = zeros(Nh,HDim-1);
CPSObj.W = zeros(VDim-1,HDim-1,Nv,Nh);
CPSObj.Theta = zeros(Nh,HDim-1);

CPSObj.ati = zeros(1,VDim-1);
CPSObj.bti = zeros(Alpha,HDim-1);
CPSObj.Wti = zeros(VDim-1,HDim-1,Nv,Alpha);

for v = 1:(VDim-1)
    CPSObj.ati(v) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h=1:(Alpha * (HDim-1))
    CPSObj.bti(h) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for k = 1:(Alpha * Nv * (VDim-1)*(HDim-1))
    CPSObj.Wti(k) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end

% Repackage the ati, bti and Wti in the necessary CPS form.
for v = 1:Nv
    CPSObj.a(v,:) = CPSObj.ati;
end
for a = 1:Alpha
    for b = 1:numel(BondMap)
        HInd = b + (a-1)*Ntr; CPSObj.b(HInd,:) = CPSObj.bti(a,:);
        for v = 1:Nv
            VInd = BondMap{b}(v);
            CPSObj.W(:,:,VInd,HInd) = CPSObj.Wti(:,:,v,a);
        end
    end
end

CPSObj.OptInds = ones(CPSObj.Np,1);
