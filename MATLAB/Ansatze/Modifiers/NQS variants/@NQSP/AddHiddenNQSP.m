% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSP(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, W and Theta.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
% interactions. Max value VDim-1.
% - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
% interactions. Max value HDim-1.
% - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
% HOrder) + (Nv x VOrder)(Nh x HOrder)
% - NQS.a = (Nv x VOrder) matrix - visible site biases.
% - NQS.b = (Nh x HOrder) matrix - hidden site bias.
% - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.VisVec = (Nv x 1) vector - visible occupancies.
% ---------------------------------

% Params requires field NhP, b, W, nphs, nmag.

NhP = round(Params.NhP); % Require integer NhP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh0 = NQSObj.Nh; % Number of "hidden" spins.
VOrder = NQSObj.VOrder; HOrder = NQSObj.HOrder;

b = NQSObj.b; W = NQSObj.W;

if NhP < 0
    if abs(NhP) >= Nh0
        error('Proposed action removes all hidden units from NQS object.');
    else
        NhF = Nh0 + NhP;
        bF = b(1:NhF,:); WF = W(1:NhF,:,:,:);
    end
else
    NhF = Nh0 + NhP;
    bF = zeros(NhP,HOrder); WF = zeros(NhP,Nv,HOrder,VOrder);
    for p = 1:NhP
        for ho = 1:HOrder
            bF(p,ho) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            for v = 1:Nv
                for vo = 1:VOrder
                    WF(p,v,ho,vo) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
                end
            end
        end
    end
    bF = [b; bF]; WF = [W; WF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = Nv + NhF + NhF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.b = bF; NQSObj.W = WF; NQSObj.OptInds = ones(NQSObj.Np,1);
end