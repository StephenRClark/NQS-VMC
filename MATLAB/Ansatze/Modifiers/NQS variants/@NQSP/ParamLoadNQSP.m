% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQSP(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
% interactions. Max value VDim-1.
% - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
% interactions. Max value HDim-1.
% - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
% HOrder) + (Nv x VOrder)(Nh x HOrder)
% - NQS.a = (Nv x VOrder) matrix - visible site biases.
% - NQS.b = (Nh x HOrder) matrix - hidden site bias.
% - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.VisVec = (Nv x 1) vector - visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x VOrder) x 1 for d/da.
% - (Nh x HOrder) x 1 for d/db.
% - (Nh x Nv) x (HOrder x VOrder) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VOrder = NQSObj.VOrder; % Number of "visible" spins.
Nh = NQSObj.Nh; HOrder = NQSObj.HOrder; % Number of "hidden" spins.

% Unpack the changes in parameters of the NQS:
da = P(1:Nv*VOrder);
db = P((1:Nh*HOrder)+Nv*VOrder);
dW = reshape(P((1:Nh*Nv*HOrder*VOrder)+Nv*VOrder+Nh*HOrder),Nh,Nv,HOrder,VOrder);

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

NQSObj.OptInds = [(real(P)~=0), (imag(P)~=0)]; % Assume the non-zero parameters are intended to be varied.

end