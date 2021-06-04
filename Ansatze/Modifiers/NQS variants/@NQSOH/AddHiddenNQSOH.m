% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSOH(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, W and Theta.
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

% Params requires field NhP, b, W, nphs, nmag.

NhP = round(Params.NhP); % Require integer NhP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
Nh0 = NQSObj.Nh; % Number of "hidden" spins.

b = NQSObj.b; W = NQSObj.W;

if NhP < 0
    if abs(NhP) >= Nh0
        error('Proposed action removes all hidden units from NQS object.');
    else
        NhF = Nh0 + NhP;
        bF = b(1:NhF); WF = W(1:NhF,:);
    end
else
    NhF = Nh0 + NhP;
    bF = zeros(NhP,1); WF = zeros(NhP,VDim*Nv);
    for p = 1:NhP
        bF(p) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:VDim*Nv
            WF(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    bF = [b; bF]; WF = [W; WF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = VDim*Nv + NhF + NhF*VDim*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.b = bF; NQSObj.W = WF; NQSObj.OptInds = ones(NQSObj.Np,1);

end