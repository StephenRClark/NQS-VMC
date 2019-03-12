% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQS(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% % ---------------------------------
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
da = dP(1:Nv);
db = dP((Nv+1):(Nv+Nh));
dW = reshape(dP((Nv+Nh+1):(Nv+Nh+Nv*Nh)),Nh,Nv);

% Apply updates to the ansatz:
NQSObj.a = NQSObj.a + da;
NQSObj.b = NQSObj.b + db;
NQSObj.W = NQSObj.W + dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.a(isinf(NQSObj.a)) = 0;
NQSObj.a(isnan(NQSObj.a)) = 0;
ind = abs(real(NQSObj.a))>cap;
NQSObj.a(ind) = sign(real(NQSObj.a(ind)))*cap + 1i*imag(NQSObj.a(ind));

NQSObj.b(isinf(NQSObj.b)) = 0;
NQSObj.b(isnan(NQSObj.b)) = 0;
ind = abs(real(NQSObj.b))>cap;
NQSObj.b(ind) = sign(real(NQSObj.b(ind)))*cap + 1i*imag(NQSObj.b(ind));

NQSObj.W(isinf(NQSObj.W)) = 0;
NQSObj.W(isnan(NQSObj.W)) = 0;
ind = abs(real(NQSObj.W))>cap;
NQSObj.W(ind) = sign(real(NQSObj.W(ind)))*cap + 1i*imag(NQSObj.W(ind));
