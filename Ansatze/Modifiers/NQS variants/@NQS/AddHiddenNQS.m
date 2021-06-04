% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQS(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, W and Theta.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------

% Params requires field NhP, b, W, nphs, nmag.

NhP = round(Params.NhP); % Require integer NhP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
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
    bF = zeros(NhP,1); WF = zeros(NhP,Nv);
    for p = 1:NhP
        bF(p) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WF(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    bF = [b; bF]; WF = [W; WF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = Nv + NhF + NhF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.b = bF; NQSObj.W = WF; NQSObj.OptInds = ones(NQSObj.Np,1);