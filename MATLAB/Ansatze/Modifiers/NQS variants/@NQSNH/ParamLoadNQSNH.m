% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQSNH(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + Nv*Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.B = (Nh x 1) vector - hidden site square bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/db.
% - (Nh x 1) for d/dB.
% - (Nh*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

% Unpack the changes in parameters of the NQS:
da = P(1:Nv);
dA = P((1:Nv) + Nv);
db = P((1:Nh) + 2*Nv);
dB = P((1:Nh) + 2*Nv + Nh);
dW = reshape(P((1:(Nv*Nh)) + 2*(Nv + Nh)),Nh,Nv);

% Apply updates to the ansatz:
NQSObj.a = NQSObj.a + da;
NQSObj.A = NQSObj.A + dA;
NQSObj.b = NQSObj.b + db;
NQSObj.B = NQSObj.B + dB;
NQSObj.W = NQSObj.W + dW;

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

NQSObj.A(isinf(NQSObj.A)) = 0;
NQSObj.A(isnan(NQSObj.A)) = 0;
ind = abs(NQSObj.A)>cap;
NQSObj.A(ind) = sign(NQSObj.A(ind))*cap;

NQSObj.B(isinf(NQSObj.B)) = 0;
NQSObj.B(isnan(NQSObj.B)) = 0;
ind = abs(NQSObj.B)>cap;
NQSObj.B(ind) = sign(NQSObj.B(ind))*cap;

NQSObj.W(isinf(NQSObj.W)) = 0;
NQSObj.W(isnan(NQSObj.W)) = 0;
ind = abs(NQSObj.W)>cap;
NQSObj.W(ind) = sign(NQSObj.W(ind))*cap;

NQSObj.OptInds = (P~=0); % Assume the non-zero parameters are intended to be varied.

end