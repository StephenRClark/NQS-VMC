% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSOHTI(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim + Alpha + (VDim*Nv * Alpha).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% Properties added with translation invariance:
% - NQS.ati = (VDim x 1) vector - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.Wv = (Alpha x VDim*Nv) matrix - reduced parameter set for TI.
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (VDim x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv*VDim x 1) for d/dW.
% Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; Alpha = NQSObj.Alpha;
GraphObj = NQSObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

dP = dP.*NQSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = dP(1:(VDim));
db = dP((1:Alpha)+VDim);
dW = reshape(dP((1:(Alpha*VDim*Nv))+Alpha+VDim),VDim*Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.ati = NQSObj.ati + da;
NQSObj.bti = NQSObj.bti + db;
NQSObj.Wv = NQSObj.Wv + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.ati(isinf(NQSObj.ati)) = 0;
NQSObj.ati(isnan(NQSObj.ati)) = 0;
ind = abs(real(NQSObj.ati))>cap;
NQSObj.ati(ind) = sign(real(NQSObj.ati(ind)))*cap + 1i*imag(NQSObj.ati(ind));
ind = abs(imag(NQSObj.ati))>pi;
NQSObj.ati(ind) = real(NQSObj.ati(ind)) + 1i*(mod(imag(NQSObj.ati(ind))+pi,2*pi)-pi);

NQSObj.bti(isinf(NQSObj.bti)) = 0;
NQSObj.bti(isnan(NQSObj.bti)) = 0;
ind = abs(real(NQSObj.bti))>cap;
NQSObj.bti(ind) = sign(real(NQSObj.bti(ind)))*cap + 1i*imag(NQSObj.bti(ind));
ind = abs(imag(NQSObj.bti))>pi;
NQSObj.bti(ind) = real(NQSObj.bti(ind)) + 1i*(mod(imag(NQSObj.bti(ind))+pi,2*pi)-pi);

NQSObj.Wv(isinf(NQSObj.Wv)) = 0;
NQSObj.Wv(isnan(NQSObj.Wv)) = 0;
ind = abs(real(NQSObj.Wv))>cap;
NQSObj.Wv(ind) = sign(real(NQSObj.Wv(ind)))*cap + 1i*imag(NQSObj.Wv(ind));
ind = abs(imag(NQSObj.Wv))>pi;
NQSObj.Wv(ind) = real(NQSObj.Wv(ind)) + 1i*(mod(imag(NQSObj.Wv(ind))+pi,2*pi)-pi);

% Repackage the parameters in the necessary form.
NQSObj.a = reshape(NQSObj.ati * ones(1,Nv),VDim*Nv,1);
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bti(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng));
                for v = 1:VDim
                    NInd = v + VDim*(n-1); WInd = v + VDim*(VInd-1);
                    % Account for enlarged lattices where Nv = Ns x Ng.
                    NQSObj.W(b+(a-1)*Ntr,WInd) = NQSObj.Wv(a,NInd);
                end
            end
        end
    end
end

end