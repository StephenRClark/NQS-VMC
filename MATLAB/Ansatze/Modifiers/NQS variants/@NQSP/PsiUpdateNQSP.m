% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSP(NQSObj,P)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSP Modifier:
% - NQSP.Nv = number of "visible" units.
% - NQSP.Nh = number of "hidden" units.
% - NQSP.Np = number of parameters in the ansatz = (Nsl x VOrder) + (Alpha x
% HOrder) + (Nv x VOrder)(Alpha x HOrder)
% - NQSP.VDim = dimension of the visible units.
% - NQSP.HDim = dimension of the hidden units.
% - NQSP.VOrder = highest power of visible unit interactions. Max value VDim-1.
% - NQSP.HOrder = highest power of hidden unit interactions. Max value HDim-1.
% - NQSP.a = (Nv x VOrder) matrix - visible site biases.
% - NQSP.av = (Nsl x VOrder) matrix - visible bias parameters
% - NQSP.b = (Nh x HOrder) matrix - hidden site bias.
% - NQSP.bv = (Alpha x HOrder) matrix - hidden bias parameters.
% - NQSP.W = (Nh x Nv x HOrder x VOrder) array - hidden-visible coupling terms.
% - NQSP.Wm = (Alpha x Nv x HOrder x VOrder) array - hidden-visible coupling parameters
% - NQSP.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSP.Theta = (Nh x HOrder) matrix - effective angles by hidden order.
% - NQSP.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSP.Rescale = flag for visible unit rescaling to [0 1] interval.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x VOrder) x 1 for d/da. Group by Sublattice > Visible order
% > [sl,vo], [sl, vo+1] ... [sl+1, vo]
% - (Alpha x HOrder) x 1 for d/db. Group by Alpha > Hidden order
% > [al, ho], [al, ho+1] ... [al+1, ho]
% - (Alpha x Nv) x (HOrder x VOrder) for d/dW. Group by Alpha > Position > Hidden order > Visible order
% > [al,v,ho,vo], [al,v,ho,vo+1] ... [al,v,ho+1,vo] ... [al,v+1,ho,vo] ... [al+1,v,ho,vo]
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VOrder = NQSObj.VOrder; % Number and dimension of "visible" units.
Alpha = NQSObj.Alpha; HOrder = NQSObj.HOrder; % Density and dimension of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

P = real(P).*NQSObj.OptInds(:,1) + 1i*imag(P).*NQSObj.OptInds(:,2); % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = reshape(P(1:Nsl*VOrder),VOrder,Nsl).';
db = reshape(P((1:Alpha*HOrder)+Nsl*VOrder),HOrder,Alpha).';
dW = permute(reshape(P((1:Alpha*Nv*HOrder*VOrder)+Nsl*VOrder+Alpha*HOrder),...
    VOrder,HOrder,Nv,Alpha),[4 3 2 1]);

% Apply updates to the ansatz:
NQSObj.av = NQSObj.av + da;
NQSObj.bv = NQSObj.bv + db;
NQSObj.Wm = NQSObj.Wm + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.av(isinf(NQSObj.av)) = 0;
NQSObj.av(isnan(NQSObj.av)) = 0;
ind = abs(NQSObj.av)>cap;
NQSObj.av(ind) = sign(NQSObj.av(ind))*cap;

NQSObj.bv(isinf(NQSObj.bv)) = 0;
NQSObj.bv(isnan(NQSObj.bv)) = 0;
ind = abs(NQSObj.bv)>cap;
NQSObj.bv(ind) = sign(NQSObj.bv(ind))*cap;

NQSObj.Wm(isinf(NQSObj.Wm)) = 0;
NQSObj.Wm(isnan(NQSObj.Wm)) = 0;
ind = abs(NQSObj.Wm)>cap;
NQSObj.Wm(ind) = sign(NQSObj.Wm(ind))*cap;

% Repackage the av, bv and Wm to usual NQS form.
for n = 1:Nv
    NQSObj.a(n,:) = NQSObj.av(SLInds(n),:);
end
% Constructing shift invariant W matrix.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr,:) = NQSObj.bv(al,:);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                for ho = 1:HOrder
                    for vo = 1:VOrder
                        % Account for enlarged lattices where Nv = Ns x Ng.
                        NQSObj.W(b+(al-1)*Ntr,VInd,vo,ho) = NQSObj.Wm(al,n,ho,vo);
                    end
                end
            end
        end
    end
end

end