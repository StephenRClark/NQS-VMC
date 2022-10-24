% --- General NQS wave function update function ---

function NQSObj = PsiUpdateNQSSX(NQSObj,dP)
% This function updates the NQS parameters of the ansatz from a vector of
% parameters P.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv x 1) for d/da.
% - (Nv x 1) for d/dA.
% - (Nh x 1) for d/db.
% - (Nh x 1) for d/dB.
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
db = dP((1:Nh) + 2*Nv);
dB = dP((1:Nh) + 2*Nv + Nh);
dW = reshape(dP((1:(Nv*Nh)) + 2*Nv + 2*Nh),Nh,Nv);
dX = reshape(dP((1:(Nv*Nh)) + 2*Nv + 2*Nh + Nv*Nh),Nh,Nv);

% Apply updates to the ansatz:
NQSObj.a = NQSObj.a + da;
NQSObj.A = NQSObj.A + dA;
NQSObj.b = NQSObj.b + db;
NQSObj.B = NQSObj.B + dB;
NQSObj.W = NQSObj.W + dW;
NQSObj.X = NQSObj.X + dX;

cap = NQSObj.ParamCap;

% Sanity check the values of the ansatz:
NQSObj.a = ParamCheck(NQSObj.a,cap);
NQSObj.b = ParamCheck(NQSObj.b,cap);
NQSObj.A = ParamCheck(NQSObj.A,cap);
NQSObj.B = ParamCheck(NQSObj.B,cap);
NQSObj.W = ParamCheck(NQSObj.W,cap);
NQSObj.X = ParamCheck(NQSObj.X,cap);

end