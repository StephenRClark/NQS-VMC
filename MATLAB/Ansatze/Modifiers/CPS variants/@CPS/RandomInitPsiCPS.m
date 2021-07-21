function [CPSObj] = RandomInitPsiCPS(CPSObj,Params)
Nv = CPSObj.Nv; % Number of "visible" spins.
Nh = CPSObj.Nh; % Number of "hidden" spins.

CPSObj.Np = 2 * Nv + 2 * Nh + 4 * (Nv * Nh); % The number of variational parameters.
CPSObj.a = zeros(Nv * 2,1); 
CPSObj.b = zeros(Nh * 2,1); 
CPSObj.W = zeros(Nv * Nh * 4,1);
CPSObj.Theta = zeros(Nh * 2,1);
for v = 1:(Nv * 2)
  CPSObj.a(v) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for h=1:(Nh * 2)
  CPSObj.b(h) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end
for k = 1:(Nv * Nh * 4)
  CPSObj.W(k) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end

CPSObj.OptInds = ones(CPSObj.Np,1);
