% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSMH(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P for the translation invariant structure.
% ---------------------------------
% Format for NQS Modifier object with multiplon-holon interactions:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nv + 2*Nh + 2*Nv.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.BH = (Nh x 1) vector - hidden holon bias.
% - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
% - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
% - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
% - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
% - NQS.Hv = (Nv x 1) vector - vector of visible holons.
% - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/dBH.
% - (Nh x 1) for d/dBM.
% - (Nh*Nv x 1) for d/dW.
% - (Nh*Nv x 1) for d/dX.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

dP = dP.*NQSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = dP(1:Nv);
dA = dP((1:Nv) + Nv);
dBH = dP((1:Nh) + 2*Nv);
dBM = dP((1:Nh) + 2*Nv + Nh);
dW = reshape(dP((1:(Nv*Nh)) + 2*(Nv + Nh)),Nh,Nv);
dX = reshape(dP((1:(Nv*Nh)) + 2*(Nv + Nh) + Nv*Nh),Nh,Nv);

% Apply updates to the ansatz:
NQSObj.a = NQSObj.a + da;
NQSObj.A = NQSObj.A + dA;
NQSObj.BH = NQSObj.BH + dBH;
NQSObj.BM = NQSObj.BM + dBM;
NQSObj.W = NQSObj.W + dW;
NQSObj.X = NQSObj.X + dX;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.a(isinf(NQSObj.a)) = 0;
NQSObj.a(isnan(NQSObj.a)) = 0;
ind = abs(NQSObj.a)>cap;
NQSObj.a(ind) = sign(NQSObj.a(ind))*cap;

NQSObj.A(isinf(NQSObj.A)) = 0;
NQSObj.A(isnan(NQSObj.A)) = 0;
ind = abs(NQSObj.A)>cap;
NQSObj.A(ind) = sign(NQSObj.A(ind))*cap;

NQSObj.BH(isinf(NQSObj.BH)) = 0;
NQSObj.BH(isnan(NQSObj.BH)) = 0;
ind = abs(NQSObj.BH)>cap;
NQSObj.BH(ind) = sign(NQSObj.BH(ind))*cap;

NQSObj.BM(isinf(NQSObj.BM)) = 0;
NQSObj.BM(isnan(NQSObj.BM)) = 0;
ind = abs(NQSObj.BM)>cap;
NQSObj.BM(ind) = sign(NQSObj.BM(ind))*cap;

NQSObj.W(isinf(NQSObj.W)) = 0;
NQSObj.W(isnan(NQSObj.W)) = 0;
ind = abs(NQSObj.W)>cap;
NQSObj.W(ind) = sign(NQSObj.W(ind))*cap;

NQSObj.X(isinf(NQSObj.X)) = 0;
NQSObj.X(isnan(NQSObj.X)) = 0;
ind = abs(NQSObj.X)>cap;
NQSObj.X(ind) = sign(NQSObj.X(ind))*cap;