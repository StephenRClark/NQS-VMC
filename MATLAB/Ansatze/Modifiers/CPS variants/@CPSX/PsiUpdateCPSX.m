% --- General CPS wave function update function ---

% - Original by X Fang, updated by M Pei.

function CPSObj = PsiUpdateCPSX(CPSObj,dP)
% This function updates the CPS parameters of the ansatz from a vector of
% parameters P.
% X indicates exponential parameterisation.
% ---------------------------------
% Format for CPSX Modifier object:
% - CPSX.Nv = number of "visible" spins.
% - CPSX.Nh = number of "hidden" spins.
% - CPSX.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPSX.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPSX.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPSX.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPSX.HDim = 3 - this version features fixed hidden unit dimension.
% - CPSX.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPSX.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPSX.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPSX.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPSX.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nv*(VDim-1) x 1) for d/da.
% Arranged ([i,1],...,[i,VDim],[i+1,1],...).
% - (Nh*(HDim-1) x 1) for d/db.
% Arranged ([i,1],...,[i,HDim],[i+1,1],...).
% - ((Nh*Nv*(VDim-1)*(HDim-1)) x 1) for d/dW.
% Arranged
% ([v(i),h(j),i,j],[v(i)+1,h(j),i,j],...,[v(i),h(j)+1,i,j],
%           ...,[v(i+1),h(j),i+1,j],...[v(i),h(j+1),i,j+1],...).
% ---------------------------------


% Make local copies to reduce notation in code below.
Nv = CPSObj.Nv; VDim = CPSObj.VDim; % Number of "visible" spins and visible dimension.
Nh = CPSObj.Nh; HDim = CPSObj.HDim; % Number of "hidden" spins and hidden dimension.

dP = dP.*CPSObj.OptInds; % Zeroes out any undesired parameter changes.

% Unpack the changes in parameters of the NQS:
da = reshape(dP(1:(Nv*(VDim-1))),Nv,VDim-1);
db = reshape(dP((1:(Nh*(HDim-1))) + Nv*(VDim-1)),Nh,HDim-1);
dW = reshape(dP((1:(Nh*Nv*(VDim-1)*(HDim-1))) + Nv*(VDim-1) + Nh*(HDim-1)),...
    VDim-1,HDim-1,Nv,Nh);

% Apply updates to the ansatz:
CPSObj.a = CPSObj.a + da;
CPSObj.b = CPSObj.b + db;
CPSObj.W = CPSObj.W + dW;

cap = CPSObj.ParamCap;

% Sanity check the values of the ansatz:
CPSObj.a(isinf(CPSObj.a)) = 0;
CPSObj.a(isnan(CPSObj.a)) = 0;
ind = abs(CPSObj.a)>cap;
CPSObj.a(ind) = sign(CPSObj.a(ind))*cap;

CPSObj.b(isinf(CPSObj.b)) = 0;
CPSObj.b(isnan(CPSObj.b)) = 0;
ind = abs(CPSObj.b)>cap;
CPSObj.b(ind) = sign(CPSObj.b(ind))*cap;

CPSObj.W(isinf(CPSObj.W)) = 0;
CPSObj.W(isnan(CPSObj.W)) = 0;
ind = abs(CPSObj.W)>cap;
CPSObj.W(ind) = sign(CPSObj.W(ind))*cap;
end
