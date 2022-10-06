% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQS(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.Alpha = number of unique coupling sets or "hidden unit density".
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.av = (Nsl x 1) vector - visible bias parameters.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% Arranged [a, v], [a, v+1] ... [a+1, v] ...
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
db = P((1:Alpha)+Nsl);
dW = reshape(P((1:(Alpha*Nv))+Nsl+Alpha),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.av = da;
NQSObj.bv = db;
NQSObj.Wm = dW;

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
NQSObj.bv(ind) = sign(real(NQSObj.bv(ind)))*cap + 1i*imag(NQSObj.bv(ind));
ind = abs(imag(NQSObj.bv))>pi;
NQSObj.bv(ind) = real(NQSObj.bv(ind)) + 1i*(mod(imag(NQSObj.bv(ind))+pi,2*pi)-pi);

NQSObj.Wm(isinf(NQSObj.Wm)) = 0;
NQSObj.Wm(isnan(NQSObj.Wm)) = 0;
ind = abs(real(NQSObj.Wm))>cap;
NQSObj.Wm(ind) = sign(real(NQSObj.Wm(ind)))*cap + 1i*imag(NQSObj.Wm(ind));
ind = abs(imag(NQSObj.Wm))>pi;
NQSObj.Wm(ind) = real(NQSObj.Wm(ind)) + 1i*(mod(imag(NQSObj.Wm(ind))+pi,2*pi)-pi);

% Repackage the ati, bti and Wv to usual NQS form.
for n = 1:Nv
    NQSObj.a(n) = NQSObj.av(SLInds(n));
end
% Constructing shift invariant W matrix.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
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

% Set optimisation indicies.
NQSObj.OptInds = [(real(P)~=0), (imag(P)~=0)]; 

end