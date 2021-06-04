% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSOH(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
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

dP = dP.*NQSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = dP(1:(VDim*Nv));
db = dP((1:Nh)+(VDim*Nv));
dW = reshape(dP((1:(Nh*VDim*Nv))+Nh+VDim*Nv),VDim*Nv,Nh).';

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
ind = abs(imag(NQSObj.a))>pi;
NQSObj.a(ind) = real(NQSObj.a(ind)) + 1i*(mod(imag(NQSObj.a(ind))+pi,2*pi)-pi);

NQSObj.b(isinf(NQSObj.b)) = 0;
NQSObj.b(isnan(NQSObj.b)) = 0;
ind = abs(real(NQSObj.b))>cap;
NQSObj.b(ind) = sign(real(NQSObj.b(ind)))*cap + 1i*imag(NQSObj.b(ind));
ind = abs(imag(NQSObj.b))>pi;
NQSObj.b(ind) = real(NQSObj.b(ind)) + 1i*(mod(imag(NQSObj.b(ind))+pi,2*pi)-pi);

NQSObj.W(isinf(NQSObj.W)) = 0;
NQSObj.W(isnan(NQSObj.W)) = 0;
ind = abs(real(NQSObj.W))>cap;
NQSObj.W(ind) = sign(real(NQSObj.W(ind)))*cap + 1i*imag(NQSObj.W(ind));
ind = abs(imag(NQSObj.W))>pi;
NQSObj.W(ind) = real(NQSObj.W(ind)) + 1i*(mod(imag(NQSObj.W(ind))+pi,2*pi)-pi);

end