% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSTI(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P for the translation invariant structure.
% ---------------------------------
% Format for NQS Modifier object with translation invariance:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nh + Alpha + 1.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% Properties added with translation invariance:
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dWv.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Alpha = NQSObj.Alpha; % Hidden unit density, needs to be integer
Ns = Nv/Ng; % Determine NQS 'lattice' size factor versus the actual Graph.

dP = dP.*NQSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = dP(1:Ns);
db = dP((1:Alpha)+Ns);
dW = reshape(dP((Alpha+Ns+1):end),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.ati = NQSObj.ati + da;
NQSObj.bti = NQSObj.bti + db;
NQSObj.Wv = NQSObj.Wv + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.ati(isinf(NQSObj.ati)) = 0;
NQSObj.ati(isnan(NQSObj.ati)) = 0;
ind = abs(NQSObj.ati)>cap;
NQSObj.ati(ind) = sign(NQSObj.ati(ind))*cap;

NQSObj.bti(isinf(NQSObj.bti)) = 0;
NQSObj.bti(isnan(NQSObj.bti)) = 0;
ind = abs(NQSObj.bti)>cap;
NQSObj.bti(ind) = sign(NQSObj.bti(ind))*cap;

NQSObj.Wv(isinf(NQSObj.Wv)) = 0;
NQSObj.Wv(isnan(NQSObj.Wv)) = 0;
ind = abs(NQSObj.Wv)>cap;
NQSObj.Wv(ind) = sign(NQSObj.Wv(ind))*cap;

% Repackage the ati, bti and Wv to usual NQS form.
for s = 1:Ns
    NQSObj.a((1:Ng)+(s-1)*Ng) = NQSObj.ati(s) * ones(Ng,1);
end

% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wv(a,n);
            end
        end
    end
end

end