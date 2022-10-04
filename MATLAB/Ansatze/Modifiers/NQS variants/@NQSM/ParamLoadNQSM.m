% --- General NQS wave function update function ---

function NQSObj = ParamLoadNQSM(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.VDim = dimensions of the visible units.
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (3*Nsl x 1) for d/da.
% Arranged [sl, H], [sl+1, H] ... [sl, D] ... [sl, M] ...
% - (Alpha x 1) for d/db.
% - (Alpha*Nv x 1) for d/dW.
% Arranged [a, v], [a, v+1] ... [a+1, v] ...
% - (Alpha*Nv x 1) for d/dX.
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
da = P(1:(3*Nsl));
db = P((1:Alpha)+(3*Nsl));
dW = reshape(P((1:(Alpha*Nv))+Alpha+3*Nsl),Nv,Alpha).';
dX = reshape(P((1:(Alpha*Nv))+Alpha*(1+Nv)+3*Nsl),Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.av = da;
NQSObj.bv = db;
NQSObj.Wm = dW;
NQSObj.Xm = dX;

cap = NQSObj.ParamCap;

NQSObj.Xm(isinf(NQSObj.Xm)) = 0;
NQSObj.Xm(isnan(NQSObj.Xm)) = 0;
ind = abs(real(NQSObj.Xm))>cap;
NQSObj.Xm(ind) = sign(real(NQSObj.Xm(ind)))*cap + 1i*imag(NQSObj.Xm(ind));
ind = abs(imag(NQSObj.Xm))>pi;
NQSObj.Xm(ind) = real(NQSObj.Xm(ind)) + 1i*(mod(imag(NQSObj.Xm(ind))+pi,2*pi)-pi);

% Repackage the parameters in the necessary form.
a_h = zeros(Nv,1); a_d = zeros(Nv,1); a_m = zeros(Nv,1);
for s = 1:Nsl
    a_h(SLInds==s) = NQSObj.av(s);
    a_d(SLInds==s) = NQSObj.av(s+Nsl);
    a_m(SLInds==s) = NQSObj.av(s+2*Nsl);
end
NQSObj.a = [a_h; a_d; a_m];
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for bd = 1:Ntr
        for n = 1:Nv
            if BondMap{bd}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{bd}(1+mod(n-1,Ng));
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(bd+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
                NQSObj.X(bd+(al-1)*Ntr,VInd) = NQSObj.Xm(al,n);
            end
        end
    end
end
% Set optimisation indicies.
NQSObj.OptInds = [(real(P)~=0), (imag(P)~=0)]; % Assume the non-zero parameters are intended to be varied.

end