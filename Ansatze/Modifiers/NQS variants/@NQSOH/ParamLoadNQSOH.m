% --- General NQS wave function update function ---

function NQSObj = ParamLoadNQSOH(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQSOH Modifier object:
% - NQSOH.Nv = number of "visible" spins.
% - NQSOH.Nh = number of "hidden" spins.
% - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
% - NQSOH.VDim = dimensions of the visible units.
% - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
% - NQSOH.b = (Nh x 1) vector - hidden site bias.
% - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
% - NQSOH.Theta = (Nh x 1) vector - effective angles.
% - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (VDim*Nv x 1) for d/da.
% Arranged [v, vd], [v, vd+1], ... , [v+1, vd], ...
% - (Nh x 1) for d/db.
% - (Nh*Nv*VDim x 1) for d/dW.
% Arranged [h, v, vd], [h, v, vd+1], ... ,[h, v+1, vd], ...
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VDim = NQSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = NQSObj.Nh; % Number of "hidden" spins.

% Unpack the changes in parameters of the NQS:
da = P(1:(VDim*Nv));
db = P((1:Nh)+(VDim*Nv));
dW = reshape(P((1:(Nh*VDim*Nv))+Nh+VDim*Nv),Nh,VDim*Nv);

% Apply updates to the ansatz:
NQSObj.a = da;
NQSObj.b = db;
NQSObj.W = dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.a(isinf(NQSObj.a)) = 0;
NQSObj.a(isnan(NQSObj.a)) = 0;
ind = abs(NQSObj.a)>cap;
NQSObj.a(ind) = sign(NQSObj.a(ind))*cap;

NQSObj.b(isinf(NQSObj.b)) = 0;
NQSObj.b(isnan(NQSObj.b)) = 0;
ind = abs(NQSObj.b)>cap;
NQSObj.b(ind) = sign(NQSObj.b(ind))*cap;

NQSObj.W(isinf(NQSObj.W)) = 0;
NQSObj.W(isnan(NQSObj.W)) = 0;
ind = abs(NQSObj.W)>cap;
NQSObj.W(ind) = sign(NQSObj.W(ind))*cap;

NQSObj.OptInds = (P~=0); % Assume the non-zero parameters are intended to be varied.

end