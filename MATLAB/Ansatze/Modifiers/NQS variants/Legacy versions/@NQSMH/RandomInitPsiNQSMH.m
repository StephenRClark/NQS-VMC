% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSMH(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
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

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

NQSObj.Np = 2*Nv + 2*Nh + 2*(Nv * Nh); % The number of variational parameters.

% Initialise the storage:
NQSObj.a = zeros(Nv,1);
NQSObj.A = zeros(Nv,1);
NQSObj.BH = zeros(Nh,1);
NQSObj.BM = zeros(Nh,1);
NQSObj.W = zeros(Nh,Nv);
NQSObj.X = zeros(Nh,Nv);
NQSObj.ThetaH = zeros(Nh,1);
NQSObj.ThetaM = zeros(Nh,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
end
if (isfield(Params,'BH') == 0) && isfield(Params,'BM')
    Params.BH = Params.BM;
elseif (isfield(Params,'BM') == 0) && isfield(Params,'BH')
    Params.BM = Params.BH;
elseif (isfield(Params,'BH') == 0) && (isfield(Params,'BM') == 0) && isfield(Params,'b')
    Params.BH = Params.b; Params.BM = Params.b;
end
if isfield(Params,'X') == 0
    Params.X = Params.W;
end

for v = 1:Nv
    NQSObj.a(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
    NQSObj.A(v) = (Params.A + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.A~=0);
end
for h=1:Nh
    NQSObj.BH(h) = (Params.BH + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.BH~=0);
    NQSObj.BM(h) = (Params.BM + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.BM~=0);
end
for h = 1:Nh
    for v = 1:Nv
        NQSObj.W(h,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
        NQSObj.X(h,v) = (Params.X + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.X~=0);
    end
end

NQSObj.OptInds = ones(NQSObj.Np,1);
end