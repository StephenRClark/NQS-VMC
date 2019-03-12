% --- General NQS spin wave function random initialisation function ---

function [Ansatz] = RandomInitPsiNQSTIRI2D_Spin(Ansatz,Params)
% This function populates random initial NQS ansatz structure for
% a spin-1/2 system. The input Ansatz is assumed to have Nv, Nh, L1/L2 
% defined already. The Params structure contains information controlling 
% the form of random elements generated.
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

% Make local copies to reduce notation in code below.
Nv = Ansatz.Nv; % Number of "visible" spins.
Nh = Ansatz.Nh; % Number of "hidden" spins.
Alpha = round(Nh/Nv); % Hidden unit density, needs to be integer.
Ansatz.Np = Nh + Alpha + 1; % The number of variational parameters.
L1 = Ansatz.Dim(1); L2 = Ansatz.Dim(2); % Required for configuration shift indexing.
if L1 == L2
    SymFac = 4;
else 
    SymFac = 2;
end
Ansatz.Sym = SymFac;
AlphaRot = SymFac * Alpha; NhRot = Nh * SymFac;
Ansatz.Nh = Nv * AlphaRot;

% Initialise the storage:
Ansatz.ati = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
Ansatz.b = zeros(NhRot,1);
Ansatz.bti = zeros(Alpha,1);
Ansatz.W = zeros(NhRot,Nv);
Ansatz.Wv = zeros(Alpha,Nv);

Ansatz.Theta = zeros(NhRot,1);

for a = 1:Alpha
    Ansatz.bti(a) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    for v = 1:Nv
        Ansatz.b(v+(a-1)*Nv) = Ansatz.bti(a);
        Ansatz.Wv(a,v) = Params.c * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    end
end
Ansatz.a = Ansatz.ati * ones(Nv,1);

% Constructing shift and rotation invariant W matrix from parameters just generated.
for k = 1:Alpha
    for dx = 0:(L1-1)
        for dy = 0:(L2-1)
            HInd = (1 + SymFac*((dy * L1) + dx) + (k-1) * SymFac * Nv);
            Ansatz.b(HInd+(0:(SymFac-1))) = Ansatz.bti(k)*ones(SymFac,1);
            for m = 1:Nv
                I = mod(m + dx - 1,L1); J = mod((ceil(m/L1) + dy - 1),L2);
                VInd0 = 1 + I + L1 * J;
                VInd180 = Nv - I - L1 * J;
                Ansatz.W(HInd,VInd0) = Ansatz.Wv(k,m);
                if SymFac == 4
                    VInd90 = 1 + J + (L2-1-I)*L1;
                    VInd270 = (I+1)*L2 - J;
                    Ansatz.W(HInd+1,VInd90) = Ansatz.Wv(k,m);
                    Ansatz.W(HInd+2,VInd180) = Ansatz.Wv(k,m);
                    Ansatz.W(HInd+3,VInd270) = Ansatz.Wv(k,m);
                else
                    Ansatz.W(HInd+1,VInd180) = Ansatz.Wv(k,m);
                end
            end
        end
    end
end

Ansatz.OptInds = ones(Ansatz.Np,1);