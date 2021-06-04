% --- General CPS wave function parameter overwrite function ---

function CPSObj = ParamLoadCPSTI(CPSObj,P)
% This function replaces the CPS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for CPS Modifier object:
% - CPS.Nv = number of "visible" spins.
% - CPS.Nh = number of "hidden" spins.
% - CPS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPS.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPS.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPS.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPS.HDim = 3 - this version features fixed hidden unit dimension.
% - CPS.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPS.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPS.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPS.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPS.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% Properties added with translation invariance:
% - CPS.Alpha = number of distinct hidden unit sets.
% - CPS.ati = (1 x (VDim-1)) vector - reduced parameter set for TI.
% - CPS.bti = (Alpha x (HDim-1)) matrix - reduced parameter set for TI.
% - CPS.Wti = ((VDim-1) x (HDim-1) x Nv x Alpha) array - reduced parameter set for TI.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP and an
% updated local index vector VisInds.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - ((VDim-1) x 1) for d/da.
% Arranged ([i,1],...,[i,VDim],[i+1,1],...).
% - (Alpha*(HDim-1) x 1) for d/db.
% Arranged ([i,1],...,[i,HDim],[i+1,1],...).
% - ((Alpha*Nv*(VDim-1)*(HDim-1)) x 1) for d/dW.
% Arranged
% ([v(i),h(a),i,a],[v(i)+1,h(a),i,a],...,[v(i),h(a)+1,i,a],
%           ...,[v(i+1),h(a),i+1,a],...[v(i),h(a+1),i,a+1],...).
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = CPSObj.Nv; VDim = CPSObj.VDim; % Number of "visible" spins and visible dimension.
Alpha = CPSObj.Alpha; HDim = CPSObj.HDim; % Number of hidden stacks and hidden dimension.
BondMap = CPSObj.Graph.BondMap; Ntr = numel(BondMap);

P = P.*CPSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = reshape(P(1:(VDim-1)),1,VDim-1);
db = reshape(P((1:(Alpha*(HDim-1))) + (VDim-1)),Alpha,HDim-1);
dW = reshape(P((1:(Alpha*Nv*(VDim-1)*(HDim-1))) + (VDim-1) + Alpha*(HDim-1)),...
    VDim-1,HDim-1,Nv,Alpha);

% Apply updates to the ansatz:
CPSObj.ati = da;
CPSObj.bti = db;
CPSObj.Wti = dW;

cap = CPSObj.ParamCap;

% Sanity check the values of the ansatz:
CPSObj.ati(isinf(CPSObj.ati)) = 0;
CPSObj.ati(isnan(CPSObj.ati)) = 0;
ind = abs(CPSObj.ati)>cap;
CPSObj.ati(ind) = sign(CPSObj.ati(ind))*cap;

CPSObj.bti(isinf(CPSObj.bti)) = 0;
CPSObj.bti(isnan(CPSObj.bti)) = 0;
ind = abs(CPSObj.bti)>cap;
CPSObj.bti(ind) = sign(CPSObj.bti(ind))*cap;

CPSObj.Wti(isinf(CPSObj.Wti)) = 0;
CPSObj.Wti(isnan(CPSObj.Wti)) = 0;
ind = abs(CPSObj.Wti)>cap;
CPSObj.Wti(ind) = sign(CPSObj.Wti(ind))*cap;

% Repackage the ati, bti and Wti in the necessary CPS form.
for v = 1:Nv
    CPSObj.a(v,:) = CPSObj.ati;
end
for a = 1:Alpha
    for b = 1:numel(BondMap)
        HInd = b + (a-1)*Ntr; CPSObj.b(HInd,:) = CPSObj.bti(a,:);
        for v = 1:Nv
            VInd = BondMap{b}(v);
            CPSObj.W(:,:,VInd,HInd) = CPSObj.Wti(:,:,v,a);
        end
    end
end

CPSObj.OptInds = (P~=0); % Assume the non-zero parameters are intended to be varied.

end