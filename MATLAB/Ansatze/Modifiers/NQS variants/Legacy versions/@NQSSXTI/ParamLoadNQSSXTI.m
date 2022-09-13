% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQSSXTI(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
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
% Properties added with translation invariance:
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------
% Format for Params vector is a vertically concatenated stack:
% - (1 x 1) for a.
% - (1 x 1) for A.
% - (Alpha x 1) for b.
% - (Alpha x 1) for B.
% - (Alpha*Nv x 1) for W.
% - (Alpha*Nv x 1) for X.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Alpha = NQSObj.Alpha; % Hidden unit density, needs to be integer

% Unpack the changes in parameters of the NQS:
da = dP(1);
dA = dP(2);
db = dP((1:Alpha)+2);
dB = dP((1:Alpha)+2+Alpha);
dW = reshape(dP((1:(Nv*Alpha))+2+2*Alpha),Nv,Alpha).';
dX = reshape(dP((1:(Nv*Alpha))+2+Alpha*(2+Nv)),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.ati = da;
NQSObj.Ati = dA;
NQSObj.bti = db;
NQSObj.Bti = dB;
NQSObj.Xv = dX;
NQSObj.Wv = dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.ati = ParamCheck(NQSObj.ati,cap);
NQSObj.Ati = ParamCheck(NQSObj.Ati,cap);
NQSObj.bti = ParamCheck(NQSObj.bti,cap);
NQSObj.Bti = ParamCheck(NQSObj.Bti,cap);
NQSObj.Wv = ParamCheck(NQSObj.Wv,cap);
NQSObj.Xv = ParamCheck(NQSObj.Xv,cap);

% Repackage the ati, bti and Wv to usual NQS form.
NQSObj.a = NQSObj.ati * ones(Nv,1);
NQSObj.A = NQSObj.Ati * ones(Nv,1);

% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    NQSObj.B((1:Ntr)+(a-1)*Ntr) = NQSObj.Bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wv(a,n);
                NQSObj.X(b+(a-1)*Ntr,VInd) = NQSObj.Xv(a,n);
            end
        end
    end
end

NQSObj.OptInds = (P~=0); % Assume the non-zero parameters are intended to be varied.

end