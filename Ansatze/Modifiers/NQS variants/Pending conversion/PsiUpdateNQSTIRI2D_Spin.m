% --- General NQS spin wave function update function ---

function Ansatz = PsiUpdateNQSTIRI2D_Spin(Ansatz,P)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P for the translation and rotation invariant structure.
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
SymFac = Ansatz.Sym;
Nv = Ansatz.Nv; % Number of "visible" spins.
Nh = Ansatz.Nh/SymFac; % Number of "hidden" spins, without SymFac enhancement.
Alpha = round(Nh/Nv); % Hidden unit density without SymFac enhancement.
L1 = Ansatz.Dim(1); L2 = Ansatz.Dim(2); % Required for configuration shift indexing.

% Unpack the changes in parameters of the NQS:
da = P(1);
db = P(2:(Alpha+1));
dW = reshape(P((Alpha+2):end),Nv,Alpha).';

% Apply updates to the ansatz:
Ansatz.ati = Ansatz.ati + da;
Ansatz.bti = Ansatz.bti + db;
Ansatz.Wv = Ansatz.Wv + dW;

cap = Ansatz.ParamCap;

% Sanity check the values of the ansatz:
Ansatz.ati(isinf(Ansatz.ati)) = 0;
Ansatz.ati(isnan(Ansatz.ati)) = 0;
ind = abs(real(Ansatz.ati))>cap;
Ansatz.ati(ind) = sign(real(Ansatz.ati(ind)))*cap + 1i*imag(Ansatz.ati(ind));

Ansatz.bti(isinf(Ansatz.bti)) = 0;
Ansatz.bti(isnan(Ansatz.bti)) = 0;
ind = abs(real(Ansatz.bti))>cap;
Ansatz.bti(ind) = sign(real(Ansatz.bti(ind)))*cap + 1i*imag(Ansatz.bti(ind));

Ansatz.Wv(isinf(Ansatz.Wv)) = 0;
Ansatz.Wv(isnan(Ansatz.Wv)) = 0;
ind = abs(real(Ansatz.Wv))>cap;
Ansatz.Wv(ind) = sign(real(Ansatz.Wv(ind)))*cap + 1i*imag(Ansatz.Wv(ind));

% Repackage the ati, bti and Wv to usual NQS form

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