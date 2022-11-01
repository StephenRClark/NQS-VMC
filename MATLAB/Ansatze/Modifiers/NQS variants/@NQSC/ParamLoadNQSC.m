% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQSC(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSC Modifier:
% - NQSC.Nv = number of "visible" units.
% - NQSC.Nh = number of "hidden" units.
% - NQSC.Np = number of parameters in the ansatz = 2*Alpha + 2*Alpha*Nv + 2*Nsl.
% - NQSC.a = (Nv x 1) vector - visible site bias.
% - NQSC.av = (Nsl x 1) vector - visible bias parameters.
% - NQSC.A = (Nv x 1) vector - visible site square bias.
% - NQSC.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSC.b = (Nh x 1) vector - hidden site bias.
% - NQSC.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSC.B = (Nh x 1) vector- hidden site square bias.
% - NQSC.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSC.w = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSC.wm = (Alpha x Nv) matrix - coupling parameters
% - NQSC.W = (Nh x Nv) matrix - hidden-square-visible coupling terms.
% - NQSC.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSC.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSC.HDim = dimension of the hidden units.
% - NQSC.Theta = (Nh x 1) vector - effective angles.
% - NQSC.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB
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

% Unpack the changes in parameters of the NQS:
da = P(1:Nsl);
dA = P((1:Nsl) + Nsl);
db = P((1:Alpha) + 2*Nsl);
dB = P((1:Alpha) + 2*Nsl + Alpha);
dw = reshape(P((1:(Nv*Alpha)) + 2*Nsl + 2*Alpha),Nv,Alpha).';
dW = reshape(P((1:(Nv*Alpha)) + 2*Nsl + 2*Alpha + Alpha*Nv),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.av = da;
NQSObj.Av = dA;
NQSObj.bv = db;
NQSObj.Bv = dB;
NQSObj.wm = wm;
NQSObj.Wm = dW;

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

NQSObj.wm(isinf(NQSObj.wm)) = 0;
NQSObj.wm(isnan(NQSObj.wm)) = 0;
ind = abs(real(NQSObj.wm))>cap;
NQSObj.wm(ind) = sign(real(NQSObj.wm(ind)))*cap + 1i*imag(NQSObj.wm(ind));
ind = abs(imag(NQSObj.wm))>pi;
NQSObj.wm(ind) = real(NQSObj.wm(ind)) + 1i*(mod(imag(NQSObj.wm(ind))+pi,2*pi)-pi);

NQSObj.Wm(isinf(NQSObj.Wm)) = 0;
NQSObj.Wm(isnan(NQSObj.Wm)) = 0;
ind = abs(real(NQSObj.Wm))>cap;
NQSObj.Wm(ind) = sign(real(NQSObj.Wm(ind)))*cap + 1i*imag(NQSObj.Wm(ind));
ind = abs(imag(NQSObj.Wm))>pi;
NQSObj.Wm(ind) = real(NQSObj.Wm(ind)) + 1i*(mod(imag(NQSObj.Wm(ind))+pi,2*pi)-pi);

NQSObj.OptInds = [(real(P)~=0), (imag(P)~=0)]; % Assume the non-zero parameters are intended to be varied.

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
                NQSObj.w(b+(al-1)*Ntr,VInd) = NQSObj.wm(al,n);
                NQSObj.W(b+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
            end
        end
    end
end

end