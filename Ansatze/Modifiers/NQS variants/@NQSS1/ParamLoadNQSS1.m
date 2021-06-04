% --- General NQS wave function parameter overwrite function ---

function NQSObj = ParamLoadNQSS1(NQSObj,P)
% This function replaces the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.A = (Nv x 1) vector - visible site square bias.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/db.
% - (Nh*Nv x 1) for d/dw.
% - (Nh*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

% Unpack the changes in parameters of the NQS:
da = P(1:Nv);
dA = P((1:Nv) + Nv);
db = P((1:Nh) + 2*Nv);
dw = reshape(P((1:(Nv*Nh)) + 2*Nv + Nh),Nh,Nv);
dW = reshape(P((1:(Nv*Nh)) + 2*Nv + Nh + Nv*Nh),Nh,Nv);

% Apply updates to the ansatz:
NQSObj.a = da;
NQSObj.A = dA;
NQSObj.b = db;
NQSObj.w = dw;
NQSObj.W = dW;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.a(isinf(NQSObj.a)) = 0;
NQSObj.a(isnan(NQSObj.a)) = 0;
ind = abs(real(NQSObj.a))>cap;
NQSObj.a(ind) = sign(real(NQSObj.a(ind)))*cap + 1i*imag(NQSObj.a(ind));
ind = abs(imag(NQSObj.a))>pi;
NQSObj.a(ind) = real(NQSObj.a(ind)) + 1i*(mod(imag(NQSObj.a(ind))+pi,2*pi)-pi);

NQSObj.b(isinf(NQSObj.b)) = 0;
NQSObj.b(isnan(NQSObj.b)) = 0;
ind = abs(real(NQSObj.b))>cap;
NQSObj.b(ind) = sign(real(NQSObj.b(ind)))*cap + 1i*imag(NQSObj.b(ind));
ind = abs(imag(NQSObj.b))>pi;
NQSObj.b(ind) = real(NQSObj.b(ind)) + 1i*(mod(imag(NQSObj.b(ind))+pi,2*pi)-pi);

NQSObj.A(isinf(NQSObj.A)) = 0;
NQSObj.A(isnan(NQSObj.A)) = 0;
ind = abs(real(NQSObj.A))>cap;
NQSObj.A(ind) = sign(real(NQSObj.A(ind)))*cap + 1i*imag(NQSObj.A(ind));
ind = abs(imag(NQSObj.A))>pi;
NQSObj.A(ind) = real(NQSObj.A(ind)) + 1i*(mod(imag(NQSObj.A(ind))+pi,2*pi)-pi);

NQSObj.w(isinf(NQSObj.w)) = 0;
NQSObj.w(isnan(NQSObj.w)) = 0;
ind = abs(real(NQSObj.w))>cap;
NQSObj.w(ind) = sign(real(NQSObj.w(ind)))*cap + 1i*imag(NQSObj.w(ind));
ind = abs(imag(NQSObj.w))>pi;
NQSObj.w(ind) = real(NQSObj.w(ind)) + 1i*(mod(imag(NQSObj.w(ind))+pi,2*pi)-pi);

NQSObj.W(isinf(NQSObj.W)) = 0;
NQSObj.W(isnan(NQSObj.W)) = 0;
ind = abs(real(NQSObj.W))>cap;
NQSObj.W(ind) = sign(real(NQSObj.W(ind)))*cap + 1i*imag(NQSObj.W(ind));
ind = abs(imag(NQSObj.W))>pi;
NQSObj.W(ind) = real(NQSObj.W(ind)) + 1i*(mod(imag(NQSObj.W(ind))+pi,2*pi)-pi);

NQSObj.OptInds = (P~=0); % Assume the non-zero parameters are intended to be varied.

end