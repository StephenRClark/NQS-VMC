% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSU(NQSObj,P)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" units.
% - NQSU.Nh = number of "hidden" units.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nmax*Nsl x 1) for d/da.
% Arranged [sl, vd], [sl, vd+1], ... , [sl+1, vd], ...
% - (Alpha x 1) for d/db.
% - (Alpha*Nv*Nmax x 1) for d/dW.
% Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; Nmax = NQSObj.VDim-1; % Number and dimension of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

P = real(P).*NQSObj.OptInds(:,1) + 1i*imag(P).*NQSObj.OptInds(:,2); % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = P(1:(Nmax*Nsl));
db = P((1:Alpha)+(Nmax*Nsl));
dw = reshape(P((1:(Alpha*Nmax*Nv))+Alpha+Nmax*Nsl),Nmax*Nv,Alpha).';

% Apply updates to the ansatz:
NQSObj.av = NQSObj.av + da;
NQSObj.bv = NQSObj.bv + db;
NQSObj.Wm = NQSObj.Wm + dw;

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
ind = abs(imag(NQSObj.b))>pi;
NQSObj.bv(ind) = real(NQSObj.bv(ind)) + 1i*(mod(imag(NQSObj.bv(ind))+pi,2*pi)-pi);

NQSObj.Wm(isinf(NQSObj.Wm)) = 0;
NQSObj.Wm(isnan(NQSObj.Wm)) = 0;
ind = abs(real(NQSObj.Wm))>cap;
NQSObj.Wm(ind) = sign(real(NQSObj.Wm(ind)))*cap + 1i*imag(NQSObj.Wm(ind));
ind = abs(imag(NQSObj.Wm))>pi;
NQSObj.Wm(ind) = real(NQSObj.Wm(ind)) + 1i*(mod(imag(NQSObj.Wm(ind))+pi,2*pi)-pi);

% Repackage the parameters in the necessary form.
for n = 1:Nv
    for v = 1:Nmax
        NQSObj.a(v + (n-1)*Nmax) = NQSObj.av(v+(SLInds(v)-1)*Nmax);
    end
end
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for bd = 1:Ntr
        for n = 1:Nv
            if BondMap{bd}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{bd}(1+mod(n-1,Ng));
                for v = 1:Nmax
                    NInd = v + Nmax*(n-1); WInd = v + Nmax*(VInd-1);
                    % Account for enlarged lattices where Nv = Ns x Ng.
                    NQSObj.W(bd+(al-1)*Ntr,WInd) = NQSObj.Wm(al,NInd);
                end
            end
        end
    end
end

end