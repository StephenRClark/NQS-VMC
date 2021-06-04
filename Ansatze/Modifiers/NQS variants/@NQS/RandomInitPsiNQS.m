% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQS(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.b = (1 x Nh) vector - hidden site bias.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Theta = (1 x Nh) vector - effective angles.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
Nh = NQSObj.Nh; % Number of "hidden" spins.

NQSObj.Np = Nv + Nh + (Nv * Nh); % The number of variational parameters.

% Initialise the storage:
NQSObj.a = zeros(Nv,1); 
NQSObj.b = zeros(Nh,1); 
NQSObj.W = zeros(Nh,Nv);
NQSObj.Theta = zeros(Nh,1);

for v = 1:Nv
  NQSObj.a(v) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h=1:Nh
  NQSObj.b(h) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h = 1:Nh
  for v = 1:Nv
    NQSObj.W(h,v) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
  end
end

NQSObj.OptInds = ones(NQSObj.Np,1);
