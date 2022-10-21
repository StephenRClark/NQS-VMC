% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSS1(NQSObj,P)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQSS1.Nv = number of "visible" spins.
% - NQSS1.Nh = number of "hidden" spins.
% - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
% - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSS1.a = (Nv x 1) vector - visible site bias.
% - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
% - NQSS1.A = (Nv x 1) vector - visible site square bias.
% - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSS1.b = (Nh x 1) vector - hidden site bias.
% - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
% - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
% - NQSS1.Theta = (Nh x 1) vector - effective angles.
% - NQSS1.VisVec = (Nv x 1) vector - visible occupancies vector.
% - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dw.
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

P = real(P).*NQSObj.OptInds(:,1) + 1i*imag(P).*NQSObj.OptInds(:,2); % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = P(1:Nsl);
dA = P((1:Nsl) + Nsl);
db = P((1:Alpha) + 2*Nsl);
dw = reshape(P((1:(Nv*Alpha)) + 2*Nsl + Alpha),Nv,Alpha).';
dW = reshape(P((1:(Nv*Alpha)) + 2*Nsl + Alpha*(1+Nv)),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.av = NQSObj.av + da;
NQSObj.Av = NQSObj.Av + dA;
NQSObj.bv = NQSObj.bv + db;
NQSObj.wm = NQSObj.wm + dw;
NQSObj.Wm = NQSObj.Wm + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.av(isinf(NQSObj.av)) = 0;
NQSObj.av(isnan(NQSObj.av)) = 0;
ind = abs(real(NQSObj.av))>cap;
NQSObj.av(ind) = sign(real(NQSObj.av(ind)))*cap + 1i*imag(NQSObj.av(ind));
ind = abs(imag(NQSObj.av))>pi;
NQSObj.av(ind) = real(NQSObj.av(ind)) + 1i*(mod(imag(NQSObj.av(ind))+pi,2*pi)-pi);

NQSObj.bv(isinf(NQSObj.bv)) = 0;
NQSObj.bv(isnan(NQSObj.bv)) = 0;
ind = abs(real(NQSObj.bv))>cap;
NQSObj.bv(ind) = sign(real(NQSObj.b(ind)))*cap + 1i*imag(NQSObj.bv(ind));
ind = abs(imag(NQSObj.bv))>pi;
NQSObj.bv(ind) = real(NQSObj.bv(ind)) + 1i*(mod(imag(NQSObj.bv(ind))+pi,2*pi)-pi);

NQSObj.Av(isinf(NQSObj.Av)) = 0;
NQSObj.Av(isnan(NQSObj.Av)) = 0;
ind = abs(real(NQSObj.Av))>cap;
NQSObj.Av(ind) = sign(real(NQSObj.Av(ind)))*cap + 1i*imag(NQSObj.Av(ind));
ind = abs(imag(NQSObj.A))>pi;
NQSObj.Av(ind) = real(NQSObj.Av(ind)) + 1i*(mod(imag(NQSObj.Av(ind))+pi,2*pi)-pi);

NQSObj.wm(isinf(NQSObj.wm)) = 0;
NQSObj.wm(isnan(NQSObj.wm)) = 0;
ind = abs(real(NQSObj.wm))>cap;
NQSObj.wm(ind) = sign(real(NQSObj.wm(ind)))*cap + 1i*imag(NQSObj.wm(ind));
ind = abs(imag(NQSObj.w))>pi;
NQSObj.wm(ind) = real(NQSObj.wm(ind)) + 1i*(mod(imag(NQSObj.wm(ind))+pi,2*pi)-pi);

NQSObj.Wm(isinf(NQSObj.Wm)) = 0;
NQSObj.Wm(isnan(NQSObj.Wm)) = 0;
ind = abs(real(NQSObj.Wm))>cap;
NQSObj.Wm(ind) = sign(real(NQSObj.Wm(ind)))*cap + 1i*imag(NQSObj.Wm(ind));
ind = abs(imag(NQSObj.W))>pi;
NQSObj.Wm(ind) = real(NQSObj.Wm(ind)) + 1i*(mod(imag(NQSObj.Wm(ind))+pi,2*pi)-pi);

% Repackage the ati, bti and Wv to usual NQS form.
for n = 1:Nv
    NQSObj.a(n) = NQSObj.av(SLInds(n));
    NQSObj.A(n) = NQSObj.Av(SLInds(n));
end
% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bv(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.                
                NQSObj.w(b+(a-1)*Ntr,VInd) = NQSObj.wm(a,n);
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wm(a,n);
            end
        end
    end
end

end