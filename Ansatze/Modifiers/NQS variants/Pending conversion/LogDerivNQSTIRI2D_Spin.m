% --- General NQS spin logarithmic derivative function ---

function dLogp = LogDerivNQSTIRI2D_Spin(Ansatz,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS Spin ansatz, for a spin-1/2
% configuration specifed by the structure Cfg.
% NB: Translation invariance here assumes Nh/Nv integer, and rotation
% invariance assumes at least invariance under 180 rotation
% ---------------------------------
% Format for NQS Ansatz with translation and rotation invariance in 2D:
% - Ansatz.Nv = number of "visible" spins (defined on input).
% - Ansatz.Nh = number of "hidden" spins (defined on input), enhanced by SymFac.
% - Ansatz.Np = number of parameters in the ansatz = (Nh + Alpha)/SymFac + 1. (computed here).
% - Ansatz.Dim = (L1, L2) - dimensions of the lattice.
% - Ansatz.Sym = (1 x 1) scalar - largest n-fold rotation of lattice.
% - Ansatz.a = (Nv x 1) vector - constructed from ati field.
% - Ansatz.ati = (1 x 1) scalar - reduced parameter set for TI in 2D.
% - Ansatz.b = (1 x Nh)) vector - constructed from bti field.
% - Ansatz.bti = (1 x Alpha) vector - reduced parameter set for TI in 2D.
% - Ansatz.W = (Nh x Nv) matrix - constructed from Wv field.
% - Ansatz.Wv = (Alpha x Nv) matrix - reduced parameter set for TI in 2D.
% - Ansatz.Theta = (1 x Nh) vector - NOT separated into Alpha sectors.
% - NB: NQSTI structure can be made compatible with general NQS expressions
% provided the full sets of a,b,W are reconstructed, though that will add
% some small redundancies in the information stored
% ---------------------------------
% Hidden units ordered by translation then rotated copies i.e.
% (... Translate(M) Translate(M)+Rotation(1) ... Translate(M)+Rotation(N) Translate(M+1) ...)
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (Alpha/SymFac x 1) for d/db.
% - (Nh/SymFac x 1) for d/dWv.
% ---------------------------------

% Make local copies to reduce notation in code below.
L1 = Ansatz.Dim(1); L2 = Ansatz.Dim(2); % Required for configuration shift indexing.
SymFac = Ansatz.Sym;
Nv = Ansatz.Nv; % Number of "visible" spins.
NhAct = Ansatz.Nh/SymFac; % Number of "hidden" spins without SymFac enhancement.
AlphaAct = round(NhAct/Nv); % Hidden unit density without SymFac enhancement.

Cfg_vec = FullSpinCfg(Cfg); % Build the spin configuration vector.

dLogp = zeros(NhAct + AlphaAct + 1,1); % Initialise full vector of derivatives.

if Ansatz.OptInds(1) == 1
    dLogp(1) = sum(Cfg_vec); % Insert d/da.
end

% Accounting for shift structure of W matrix requires either construction
% of shifted Theta matrix or shifted Cfg vector - the latter is done here
for a=1:AlphaAct % Derivatives need to be computed by Alpha sector
    if Ansatz.OptInds(1+a) == 1
        dLogp(1+a) = sum(tanh(Ansatz.Theta((1:SymFac*Nv) + (a-1)*SymFac*Nv))); % Insert d/db
    end
    for v = 1:Nv
        PInd = v + (a-1)*Nv;
        if Ansatz.OptInds(1+AlphaAct+PInd) == 1
            for dx = 0:(L1-1)
                for dy = 0:(L2-1)
                    I = mod(v+dx-1,L1); J = mod(ceil(v/L1) + dy - 1,L2);
                    VInd0 = 1 + I + L1*J;
                    TInd0 = 1 + SymFac*((dy*L1) + dx) + (a-1)*SymFac*Nv;
                    VInd180 = Nv - I - (J*L1);
                    if SymFac == 4
                        VInd90 = 1 + J + (L2-1-I)*L1;
                        VInd270 = (I+1)*L2 - J;
                        TInd90 = TInd0 + 1; TInd180 = TInd0 + 2; TInd270 = TInd0 + 3;
                        VInds = [VInd0; VInd90; VInd180; VInd270];
                        TInds = [TInd0; TInd90; TInd180; TInd270];
                    else
                        TInd180 = TInd0 + 1;
                        VInds = [VInd0; VInd180];
                        TInds = [TInd0; TInd180];
                    end
                    dLogp(1+AlphaAct+PInd) = dLogp(1+AlphaAct+PInd) + ( sum(Cfg_vec(VInds) .* tanh(Ansatz.Theta(TInds))) ); %Insert d/dW
                end
            end
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
