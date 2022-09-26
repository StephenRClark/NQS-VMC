% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSMHAI(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P for the translation invariant structure.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.WH = (Nh x Nv) matrix - hidden-visible HH coupling terms.
% - NQS.WM = (Nh x Nv) matrix - hidden-visible MM coupling terms.
% - NQS.XH = (Nh x Nv) matrix - hidden-visible HM coupling terms.
% - NQS.XM = (Nh x Nv) matrix - hidden-visible MH coupling terms.
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
% - NQS.WHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.WMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XHv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.XMv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (1 x 1) for d/da.
% - (1 x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB.
% - (Alpha*Nv x 1) for d/dWH.
% - (Alpha*Nv x 1) for d/dWM.
% - (Alpha*Nv x 1) for d/dXH.
% - (Alpha*Nv x 1) for d/dXM.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
Alpha = NQSObj.Alpha; % Hidden unit density, needs to be integer

dP = dP.*NQSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = dP(1);
dA = dP(2);
dBH = dP((1:Alpha)+2);
dBM = dP((1:Alpha)+2+Alpha);
dWH = reshape(dP((1:(Nv*Alpha))+2*(Alpha+1)),Nv,Alpha).';
dWM = reshape(dP((1:(Nv*Alpha))+2*(Alpha+1)+Alpha*Nv),Nv,Alpha).';
dXH = reshape(dP((1:(Nv*Alpha))+2*(Alpha+1)+2*Alpha*Nv),Nv,Alpha).';
dXM = reshape(dP((1:(Nv*Alpha))+2*(Alpha+1)+3*Alpha*Nv),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.ati = NQSObj.ati + da;
NQSObj.Ati = NQSObj.Ati + dA;
NQSObj.BHti = NQSObj.BHti + dBH;
NQSObj.BMti = NQSObj.BMti + dBM;
NQSObj.WHv = NQSObj.WHv + dWH;
NQSObj.WMv = NQSObj.WMv + dWM;
NQSObj.XHv = NQSObj.XHv + dXH;
NQSObj.XMv = NQSObj.XMv + dXM;

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

NQSObj.WHv(isinf(NQSObj.WHv)) = 0;
NQSObj.WHv(isnan(NQSObj.WHv)) = 0;
ind = abs(NQSObj.WHv)>cap;
NQSObj.WHv(ind) = sign(NQSObj.WHv(ind))*cap;

NQSObj.WMv(isinf(NQSObj.WMv)) = 0;
NQSObj.WMv(isnan(NQSObj.WMv)) = 0;
ind = abs(NQSObj.WMv)>cap;
NQSObj.WMv(ind) = sign(NQSObj.WMv(ind))*cap;

NQSObj.XHv(isinf(NQSObj.XHv)) = 0;
NQSObj.XHv(isnan(NQSObj.XHv)) = 0;
ind = abs(NQSObj.XHv)>cap;
NQSObj.XHv(ind) = sign(NQSObj.XHv(ind))*cap;

NQSObj.XMv(isinf(NQSObj.XMv)) = 0;
NQSObj.XMv(isnan(NQSObj.XMv)) = 0;
ind = abs(NQSObj.XMv)>cap;
NQSObj.XMv(ind) = sign(NQSObj.XMv(ind))*cap;

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
                NQSObj.WH(b+(a-1)*Ntr,VInd) = NQSObj.WHv(a,n);
                NQSObj.WM(b+(a-1)*Ntr,VInd) = NQSObj.WMv(a,n);
                NQSObj.XH(b+(a-1)*Ntr,VInd) = NQSObj.XHv(a,n);
                NQSObj.XM(b+(a-1)*Ntr,VInd) = NQSObj.XMv(a,n);
            end
        end
    end
end