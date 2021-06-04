% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSS1TI(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Alpha + 2 + Alpha.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% Properties added with translation invariance:
% - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
% - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
% - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
% - NQS.wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (1 x 1) for d/da.
    % - (1 x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv x 1) for d/dw.
    % - (Alpha*Nv x 1) for d/dW.
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
dw = reshape(dP((1:(Nv*Alpha))+2+Alpha),Nv,Alpha).';
dW = reshape(dP((1:(Nv*Alpha))+2+Alpha*(1+Nv)),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.ati = NQSObj.ati + da;
NQSObj.Ati = NQSObj.Ati + dA;
NQSObj.bti = NQSObj.bti + db;
NQSObj.wv = NQSObj.wv + dw;
NQSObj.Wv = NQSObj.Wv + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.ati(isinf(NQSObj.a)) = 0;
NQSObj.ati(isnan(NQSObj.a)) = 0;
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

NQSObj.Ati(isinf(NQSObj.Ati)) = 0;
NQSObj.Ati(isnan(NQSObj.Ati)) = 0;
ind = abs(real(NQSObj.Ati))>cap;
NQSObj.Ati(ind) = sign(real(NQSObj.Ati(ind)))*cap + 1i*imag(NQSObj.Ati(ind));
ind = abs(imag(NQSObj.Ati))>pi;
NQSObj.Ati(ind) = real(NQSObj.Ati(ind)) + 1i*(mod(imag(NQSObj.Ati(ind))+pi,2*pi)-pi);

NQSObj.wv(isinf(NQSObj.wv)) = 0;
NQSObj.wv(isnan(NQSObj.wv)) = 0;
ind = abs(real(NQSObj.wv))>cap;
NQSObj.wv(ind) = sign(real(NQSObj.wv(ind)))*cap + 1i*imag(NQSObj.wv(ind));
ind = abs(imag(NQSObj.wv))>pi;
NQSObj.wv(ind) = real(NQSObj.wv(ind)) + 1i*(mod(imag(NQSObj.wv(ind))+pi,2*pi)-pi);

NQSObj.Wv(isinf(NQSObj.Wv)) = 0;
NQSObj.Wv(isnan(NQSObj.Wv)) = 0;
ind = abs(real(NQSObj.Wv))>cap;
NQSObj.Wv(ind) = sign(real(NQSObj.Wv(ind)))*cap + 1i*imag(NQSObj.Wv(ind));
ind = abs(imag(NQSObj.Wv))>pi;
NQSObj.Wv(ind) = real(NQSObj.Wv(ind)) + 1i*(mod(imag(NQSObj.Wv(ind))+pi,2*pi)-pi);

% Repackage the ati, bti and Wv to usual NQS form.
NQSObj.a = NQSObj.ati * ones(Nv,1);
NQSObj.A = NQSObj.Ati * ones(Nv,1);

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
                NQSObj.w(b+(a-1)*Ntr,VInd) = NQSObj.wv(a,n);
            end
        end
    end
end

end