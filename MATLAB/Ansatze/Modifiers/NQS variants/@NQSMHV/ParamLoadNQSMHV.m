% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQSMHV(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
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
% Properties added with translation invariance:
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.BHti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.BMti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (1 x 1) for d/dA.
% - (Alpha x 1) for d/dBH.
% - (Alpha x 1) for d/dBM.
% - (Alpha*Nv x 1) for d/dW.
% - (Alpha*Nv x 1) for d/dX.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
NNInds = ReverseBond(GraphObj.Bonds); z = size(NNInds,2);
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Alpha = NQSObj.Alpha; % Hidden unit density, needs to be integer

% Unpack the changes in parameters of the NQS:
da = P(1);
dA = P(2);
dBH = P((1:Alpha)+2);
dBM = P((1:Alpha)+2+Alpha);
dW = reshape(P((1:(Nv*Alpha))+2*(Alpha+1)),Nv,Alpha).';
dXF = reshape(P((1:(z*Alpha))+2*(Alpha+1)+Alpha*Nv),z,Alpha).';
dX = - dW; dX(:,NNInds(1,:)) = dXF;

% Apply updates to the ansatz:
NQSObj.ati = NQSObj.ati + da;
NQSObj.Ati = NQSObj.Ati + dA;
NQSObj.BHti = NQSObj.BHti + dBH;
NQSObj.BMti = NQSObj.BMti + dBM;
NQSObj.Wv = NQSObj.Wv + dW;
NQSObj.Xv = NQSObj.Xv + dX;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.ati(isinf(NQSObj.ati)) = 0;
NQSObj.ati(isnan(NQSObj.ati)) = 0;
ind = abs(NQSObj.ati)>cap;
NQSObj.ati(ind) = sign(NQSObj.ati(ind))*cap;

NQSObj.Ati(isinf(NQSObj.Ati)) = 0;
NQSObj.Ati(isnan(NQSObj.Ati)) = 0;
ind = abs(NQSObj.Ati)>cap;
NQSObj.Ati(ind) = sign(NQSObj.Ati(ind))*cap;

NQSObj.BHti(isinf(NQSObj.BHti)) = 0;
NQSObj.BHti(isnan(NQSObj.BHti)) = 0;
ind = abs(NQSObj.BHti)>cap;
NQSObj.BHti(ind) = sign(NQSObj.BHti(ind))*cap;

NQSObj.BMti(isinf(NQSObj.BMti)) = 0;
NQSObj.BMti(isnan(NQSObj.BMti)) = 0;
ind = abs(NQSObj.BMti)>cap;
NQSObj.BMti(ind) = sign(NQSObj.BMti(ind))*cap;

NQSObj.Wv(isinf(NQSObj.Wv)) = 0;
NQSObj.Wv(isnan(NQSObj.Wv)) = 0;
ind = abs(NQSObj.Wv)>cap;
NQSObj.Wv(ind) = sign(NQSObj.Wv(ind))*cap;

NQSObj.Xv(isinf(NQSObj.Xv)) = 0;
NQSObj.Xv(isnan(NQSObj.Xv)) = 0;
ind = abs(NQSObj.Xv)>cap;
NQSObj.Xv(ind) = sign(NQSObj.Xv(ind))*cap;

% Repackage the ati, bti and Wv to usual NQS form.
NQSObj.a = NQSObj.ati * ones(Nv,1);
NQSObj.A = NQSObj.Ati * ones(Nv,1);

% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.BH((1:Ntr)+(a-1)*Ntr) = NQSObj.BHti(a);
    NQSObj.BM((1:Ntr)+(a-1)*Ntr) = NQSObj.BMti(a);
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