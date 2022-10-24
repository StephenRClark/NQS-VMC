% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSNH(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, B, W and Theta.
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
% ---------------------------------

% Params requires field NhP, b (B optional), W, nphs, nmag.

NhP = round(Params.NhP); % Require integer NhP.

if isfield(Params,'B') == 0
    Params.B = Params.b;
end

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh0 = NQSObj.Nh; % Number of "hidden" spins.

b = NQSObj.b; W = NQSObj.W; B = NQSObj.B;

if NhP < 0
    if abs(NhP) >= Nh0
        error('Proposed action removes all hidden units from NQS object.');
    else
        NhF = Nh0 + NhP;
        bF = b(1:NhF); BF = B(1:NhF); WF = W(1:NhF,:);
    end
else
    NhF = Nh0 + NhP;
    bF = zeros(NhP,1); BF = zeros(NhP,1); WF = zeros(NhP,Nv);
    for p = 1:NhP
        bF(p) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BF(p) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WF(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    bF = [b; bF]; BF = [B; BF]; WF = [W; WF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = 2*Nv + 2*NhF + NhF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.b = bF; NQSObj.B = BF; NQSObj.W = WF; NQSObj.OptInds = ones(NQSObj.Np,1);