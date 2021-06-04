% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSOH(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = NQSObj.Nh; % Number of "hidden" spins.

NQSObj.Np = VDim*Nv + Nh + (VDim*Nv * Nh); % The number of variational parameters.

% Initialise the storage:
NQSObj.a = zeros(VDim*Nv,1);
NQSObj.b = zeros(Nh,1);
NQSObj.W = zeros(Nh,VDim*Nv);
NQSObj.Theta = zeros(Nh,1);

for v = 1:VDim*Nv
    NQSObj.a(v) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h=1:Nh
    NQSObj.b(h) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h = 1:Nh
    for v = 1:VDim*Nv
        NQSObj.W(h,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    end
end

NQSObj.OptInds = ones(NQSObj.Np,1);
