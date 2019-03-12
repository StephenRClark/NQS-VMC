% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSNH(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
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

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

NQSObj.Np = 2*Nv + 2*Nh + (Nv * Nh); % The number of variational parameters.

% Initialise the storage:
NQSObj.a = zeros(Nv,1);
NQSObj.A = zeros(Nv,1);
NQSObj.b = zeros(Nh,1);
NQSObj.B = zeros(Nh,1);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Theta = zeros(Nh,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
end
if isfield(Params,'B') == 0
    Params.B = Params.b;
end

for v = 1:Nv
    NQSObj.a(v) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    NQSObj.A(v) = Params.A * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h=1:Nh
    NQSObj.b(h) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    NQSObj.B(h) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h = 1:Nh
    for v = 1:Nv
        NQSObj.W(h,v) = Params.c * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    end
end

NQSObj.OptInds = ones(NQSObj.Np,1);
