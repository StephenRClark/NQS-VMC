% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSMH(NQSObj,Params)
% This function adds round(NhP/Nv) new hidden units to an existing NQSObj
% (removes if negative). This will necessitate changes in Nh, Np, BH, BM,
% W, X, ThetaH and ThetaM.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nv + 2*Nh + 2*Nv.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
% - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------

% Params requires field NhP, b (B optional), W, nphs, nmag.

NhP = round(Params.NhP); % Require integer NhP.

if (isfield(Params,'BH') == 0) && isfield(Params,'BM')
    Params.BH = Params.BM;
elseif (isfield(Params,'BM') == 0) && isfield(Params,'BH')
    Params.BM = Params.BH;
elseif (isfield(Params,'BH') == 0) && (isfield(Params,'BM') == 0) && isfield(Params,'b')
    Params.BH = Params.b; Params.BM = Params.b;
end
if isfield(Params,'X') == 0
    Params.X = Params.W;
end

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh0 = NQSObj.Nh; % Number of "hidden" spins.

BH = NQSObj.BH; BM = NQSObj.BM; W = NQSObj.W; X = NQSObj.X;

if NhP < 0
    if abs(NhP) >= Nh0
        error('Proposed action removes all hidden units from NQS object.');
    else
        NhF = Nh0 + NhP;
        BHF = BH(1:NhF); BMF = BM(1:NhF); WF = W(1:NhF,:); XF = X(1:NhF,:);
    end
else
    NhF = Nh0 + NhP;
    BHF = zeros(NhP,1); BMF = zeros(NhP,1); WF = zeros(NhP,Nv); XF = zeros(NhP,Nv);
    for p = 1:NhP
        BHF(p) = Params.BH * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        BMF(p) = Params.BM * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nv
            WF(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
            XF(p,v) = Params.X * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    BHF = [BH; BHF]; BMF = [BM; BMF]; WF = [W; WF]; XF = [X; XF];
end

% Reassign all fields affected by Nh change.
NQSObj.Np = 2*Nv + 2*NhF + 2*NhF*Nv; NQSObj.Nh = NhF;
NQSObj.ThetaH = zeros(NhF,1); NQSObj.ThetaM = zeros(NhF,1);
NQSObj.BH = BHF; NQSObj.BM = BMF; NQSObj.W = WF; NQSObj.X = XF;
NQSObj.OptInds = ones(NQSObj.Np,1);