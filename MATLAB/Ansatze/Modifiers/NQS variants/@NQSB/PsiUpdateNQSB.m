% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSB(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSB Modifier:
% - NQSB.Nv = number of "visible" spins.
% - NQSB.Nh = number of "hidden" spins.
% - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSB.a = (Nv x 1) vector - visible site bias.
% - NQSB.av = (Nsl x 1) vector - visible bias parameters.
% - NQSB.A = (Nv x 1) vector - visible site square bias.
% - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSB.b = (Nh x 1) vector - hidden site bias.
% - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSB.B = (Nh x 1) vector- hidden site square bias.
% - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSB.HDim = dimension of the hidden units.
% - NQSB.Theta = (Nh x 1) vector - effective angles.
% - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Nsl = max(SLInds); Ntr = numel(BondMap); Ng = GraphObj.N; Alpha = NQSObj.Alpha;

dP = real(dP).*NQSObj.OptInds(:,1) + 1i*imag(dP).*NQSObj.OptInds(:,2); % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = P(1:Nsl);
dA = P((1:Nsl) + Nsl);
db = P((1:Alpha) + 2*Nsl);
dB = P((1:Alpha) + 2*Nsl + Alpha);
dW = reshape(P((1:(Nv*Alpha)) + 2*Nsl + 2*Alpha),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.av = NQSObj.av + da;
NQSObj.Av = NQSObj.Av + dA;
NQSObj.bv = NQSObj.bv + db;
NQSObj.Bv = NQSObj.Bv + dB;
NQSObj.Wm = NQSObj.Wm + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.av(isinf(NQSObj.av)) = 0;
NQSObj.av(isnan(NQSObj.av)) = 0;
ind = abs(real(NQSObj.av))>cap;
NQSObj.av(ind) = sign(real(NQSObj.av(ind)))*cap + 1i*imag(NQSObj.av(ind));
ind = abs(imag(NQSObj.av))>pi;
NQSObj.av(ind) = real(NQSObj.av(ind)) + 1i*(mod(imag(NQSObj.av(ind))+pi,2*pi)-pi);

NQSObj.Av(isinf(NQSObj.Av)) = 0;
NQSObj.Av(isnan(NQSObj.Av)) = 0;
ind = abs(real(NQSObj.Av))>cap;
NQSObj.Av(ind) = sign(real(NQSObj.Av(ind)))*cap + 1i*imag(NQSObj.Av(ind));
ind = abs(imag(NQSObj.Av))>pi;
NQSObj.Av(ind) = real(NQSObj.Av(ind)) + 1i*(mod(imag(NQSObj.Av(ind))+pi,2*pi)-pi);

NQSObj.bv(isinf(NQSObj.bv)) = 0;
NQSObj.bv(isnan(NQSObj.bv)) = 0;
ind = abs(real(NQSObj.bv))>cap;
NQSObj.bv(ind) = sign(real(NQSObj.bv(ind)))*cap + 1i*imag(NQSObj.bv(ind));
ind = abs(imag(NQSObj.b))>pi;
NQSObj.bv(ind) = real(NQSObj.bv(ind)) + 1i*(mod(imag(NQSObj.bv(ind))+pi,2*pi)-pi);

NQSObj.Bv(isinf(NQSObj.Bv)) = 0;
NQSObj.Bv(isnan(NQSObj.Bv)) = 0;
ind = abs(real(NQSObj.Bv))>cap;
NQSObj.Bv(ind) = sign(real(NQSObj.Bv(ind)))*cap + 1i*imag(NQSObj.Bv(ind));
ind = abs(imag(NQSObj.B))>pi;
NQSObj.Bv(ind) = real(NQSObj.Bv(ind)) + 1i*(mod(imag(NQSObj.Bv(ind))+pi,2*pi)-pi);

NQSObj.Wm(isinf(NQSObj.Wm)) = 0;
NQSObj.Wm(isnan(NQSObj.Wm)) = 0;
ind = abs(real(NQSObj.Wm))>cap;
NQSObj.Wm(ind) = sign(real(NQSObj.Wm(ind)))*cap + 1i*imag(NQSObj.Wm(ind));
ind = abs(imag(NQSObj.Wm))>pi;
NQSObj.Wm(ind) = real(NQSObj.Wm(ind)) + 1i*(mod(imag(NQSObj.Wm(ind))+pi,2*pi)-pi);

% Repackage the ati, bti and Wv to usual NQS form.
for n = 1:Nv
    NQSObj.a(n) = NQSObj.av(SLInds(n));
    NQSObj.A(n) = NQSObj.Av(SLInds(n));
end
% Constructing shift invariant W matrix.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    NQSObj.B((1:Ntr)+(al-1)*Ntr) = NQSObj.Bv(al);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(b+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
            end
        end
    end
end

end