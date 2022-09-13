% --- General NQS wave function hidden unit addition function ---

function [NQSObj] = AddHiddenNQSU(NQSObj,Params)
% This function adds NhP new hidden units to an existing NQSObj (removes if
% negative). This will n-ecessitate changes in Nh, Np, b, W and Theta.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" spins.
% - NQSU.Nh = number of "hidden" spins.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------

% Params requires field AlphaP, b, W, nphs, nmag.
AP = round(Params.AlphaP); % Require integer AlphaP.

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; Nmax = NQSObj.VDim-1; % Number of "visible" spins and visible dimension.
A0 = NQSObj.Alpha; % Number of unique coupling sets.

GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Nsl = max(SLInds); Ntr = numel(BondMap); Ng = GraphObj.N;

bv = NQSObj.bv; Wm = NQSObj.Wm;

if AP < 0
    if abs(AP) >= A0
        error('Proposed action removes all hidden units from NQS object.');
    else
        AF = A0 + AP;
        bF = bv(1:AF); WF = Wm(1:AF,:);
    end
else
    AF = A0 + AP;
    bF = zeros(AP,1); WF = zeros(AP,Nmax*Nv);
    for p = 1:AP
        bF(p) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        for v = 1:Nmax*Nv
            WF(p,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
        end
    end
    bF = [bv; bF]; WF = [Wm; WF];
end

NhF = AF*Ntr;

% Reassign all fields affected by Nh change.
NQSObj.Np = Nmax*Nsl + AF + AF*Nmax*Nv; NQSObj.Nh = NhF; NQSObj.Theta = zeros(NhF,1);
NQSObj.bv = bF; NQSObj.Wm = WF; NQSObj.Alpha = AF;

% Repackage the parameters in the necessary form.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for bd = 1:Ntr
        for n = 1:Nv
            if BondMap{bd}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{bd}(1+mod(n-1,Ng));
                for v = 1:Nmax
                    NInd = v + Nmax*(n-1); WInd = v + Nmax*(VInd-1);
                    % Account for enlarged lattices where Nv = Ns x Ng.
                    NQSObj.W(bd+(al-1)*Ntr,WInd) = NQSObj.Wm(al,NInd);
                end
            end
        end
    end
end

NQSObj.OptInds = [(real(NQSObj.av)~=0), (imag(NQSObj.av)~=0); (real(NQSObj.bv)~=0),...
    (imag(NQSObj.bv)~=0); (real(NQSObj.Wm(:))~=0), (imag(NQSObj.Wm(:))~=0)];

end