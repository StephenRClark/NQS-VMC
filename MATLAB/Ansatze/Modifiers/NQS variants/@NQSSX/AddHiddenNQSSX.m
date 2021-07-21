% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSSX(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will necessitate changes in Nh, Np, b, B, W and Theta.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Params requires field NhP, b (B optional), W, nphs, nmag.

NhP = round(Params.NhP); % Require integer NhP.

if isfield(Params,'X') == 0
    Params.X = Params.W;
end

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh0 = NQSObj.Nh; % Number of "hidden" spins.

b = NQSObj.b; B = NQSObj.B; W = NQSObj.W; X = NQSObj.X;

if NhP < 0
    if abs(NhP) >= Nh0
        error('Proposed action removes all hidden units from NQS object.');
    else
        NhF = Nh0 + NhP;
        bF = b(1:NhF); BF = B(1:NhF); XF = X(1:NhF,:); WF = W(1:NhF,:);
    end
else
    NhF = Nh0 + NhP;
    bF = zeros(NhP,1); BF = zeros(NhP,1); XF = zeros(NhP,Nv); WF = zeros(NhP,Nv);
    for p = 1:NhP
        bF(p) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BF(p) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WF(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            XF(p,v) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    bF = [b; bF]; BF = [B; BF]; XF = [X; XF]; WF = [W; WF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = 2*Nv + 2*NhF + 2*NhF*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.b = bF; NQSObj.B = NQSObj.BF; NQSObj.W = WF; NQSObj.X = XF; NQSObj.OptInds = ones(NQSObj.Np,1);
end