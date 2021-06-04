% --- General FFNN wave function random initialisation function ---

function [FFNNObj] = RandomInitPsiFFNNTI(FFNNObj,Params)
% This function populates random initial FFNN ansatz structure. The input
% FFNN is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for FFNN Modifier object:
% - FFNN.Nv = number of "visible" spins.
% - FFNN.Nh1 = number of "hidden" spins in first layer.
% - FFNN.Nh2 = number of "hidden" spins in second layer.
% - FFNN.Np = number of parameters in the ansatz = Nv + Nh1 + Nh2 + 
% - (Nv*Nh1) + (Nh1*Nh2)
% - FFNN.b1 = (Nh1 x 1) vector - first hidden site bias.
% - FFNN.b2 = (Nh2 x 1) vector - second hidden site bias.
% - FFNN.W1 = (Nh1 x Nv) matrix - hidden1-visible coupling terms.
% - FFNN.W2 = (Nh2 x Nh1) matrix - hidden2-hidden1 coupling terms.
% - FFNN.W3 = (1 x Nh2) vector - average pooling layer matrix of ones.
% - FFNN.Theta1 = (Nh1 x 1) vector - effective angles (hidden1).
% - FFNN.Theta2 = (Nh2 x 1) vector - effective angles (hidden2).
% - FFNN.h1 = (Nh1 x 1) vector - hidden unit values (first layer).
% - FFNN.h2 = (Nh2 x 1) vector - hidden unit values (second layer).
% Properties added with translation invariance:
% - FFNN.Alpha1 = number of unique W1 parameter sets
% - FFNN.b1r = (Alpha1 x 1) vector - unique 1st layer hidden biases.
% - FFNN.W1r = (Alpha1 x Nv) matrix - unique 1st layer coupling terms.
% - FFNN.W2r = (Alpha1 x Nh2) matrix - unique 1st-2nd layer coupling terms.

% Make local copies to reduce notation in code below.
Nv = FFNNObj.Nv; % Number of "visible" spins.
Ng = numel(FFNNObj.Graph.BondMap); % Number of translates required.
Nh1 = FFNNObj.Nh1; % Number of "hidden" spins in first layer.
Nh2 = FFNNObj.Nh2; % Number of "hidden" spins in second layer.

Alpha = ceil(Nh1/Nv); Nh1 = Ng*Alpha; FFNNObj.Alpha1 = Alpha;
Np = Alpha*(1+Nv+Nh2) + Nh2; FFNNObj.Np = Np;

% Initialise the storage:
% FFNNObj.a = zeros(Nv,1);
FFNNObj.b1 = zeros(Nh1,1);
FFNNObj.b2 = zeros(Nh2,1);
FFNNObj.W1 = zeros(Nh1,Nv);
FFNNObj.W2 = zeros(Nh2,Nh1);
FFNNObj.W3 = zeros(1,Nh2);
FFNNObj.Theta1 = zeros(Nh1,1);
FFNNObj.Theta2 = zeros(Nh2,1);
% Reduced parameter sets for TI:
FFNNObj.b1r = zeros(Alpha,1);
FFNNObj.W1r = zeros(Alpha,Nv);
FFNNObj.W2r = zeros(Nh2,Alpha);

for h = 1:Alpha
  FFNNObj.b1r(h) = Params.b1 * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
  for v = 1:Nv
    FFNNObj.W1r(h,v) = Params.W1 * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
  end
end

for h = 1:Nh2
  FFNNObj.b2(h) = Params.b2 * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
  for l = 1:Alpha
    FFNNObj.W2(h,l) = Params.W2 * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
  end
end

% Repackage the parameters in the necessary form.
for a = 1:Alpha
    FFNNObj.b1((1:Ng)+(a-1)*Ng) = FFNNObj.b1r(a);
    for b = 1:Ng
        H1Ind = b + (a-1)*Ng;
        for n = 1:Nv
            VInd = FFNNObj.Graph.BondMap{b}(n);
            FFNNObj.W1(H1Ind,VInd) = FFNNObj.W1r(a,n); 
        end
    end
    for n2 = 1:Nh2
        FFNNObj.W2(n2,(1:Ng)+(a-1)*Ng);
    end
end

FFNNObj.W3 = ones(1,Nh2);
FFNNObj.OptInds = ones(Np,1);
