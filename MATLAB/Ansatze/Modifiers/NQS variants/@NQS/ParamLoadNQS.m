% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQS(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (1 x Nh) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (1 x Nh) vector - effective angles.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nh x 1) for d/db.
% - (Nh*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

% Unpack the changes in parameters of the NQS:
da = P(1:Nv);
db = P((Nv+1):(Nv+Nh));
dW = reshape(P((Nv+Nh+1):(Nv+Nh+Nv*Nh)),Nh,Nv);

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