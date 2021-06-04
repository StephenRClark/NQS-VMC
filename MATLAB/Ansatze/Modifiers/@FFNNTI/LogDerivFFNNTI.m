% --- General FFNN logarithmic derivative function ---

function dLogp = LogDerivFFNNTI(FFNNObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the NQS ansatz, for a  configuration specifed
% by the structure Cfg.
% ---------------------------------
% Format for FFNN Modifier object:
% - FFNN.Nv = number of "visible" spins.
% - FFNN.Nh1 = number of "hidden" spins in first layer.
% - FFNN.Nh2 = number of "hidden" spins in second layer.
% - FFNN.Np = number of parameters in the ansatz = Alpha1 + Alpha1*Nv +
% Alpha1*Nh2 + Nh2
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
% - FFNN.W2r = (Nh2 x Alpha1) matrix - unique 1st-2nd layer coupling terms.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Alpha1 x 1) for d/db1.
% - (Nh2 x 1) for d/db2.
% - (Alpha1*Nv x 1) for d/dW1.
% - (Nh2*Alpha1 x 1) for d/dW2.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = FFNNObj.Nv; % Number of "visible" spins.
Nh1 = FFNNObj.Nh1; % Number of "hidden" spins in first layer.
Nh2 = FFNNObj.Nh2; % Number of "hidden" spins in second layer.

Alpha1 = FFNNObj.Alpha1; % You'll need to make sure Nh1/Nv is an integer.
BondMap = FFNNObj.Graph.BondMap; % Will use this to locate all the parameters.
Ng = numel(BondMap);

Cfg_vec = FFNNObj.FullCfg(Cfg); % Build the spin configuration vector.

dLogp = zeros(FFNNObj.Np,1); % Initialise full vector of derivatives.

db1 = FFNNObj.W3 * (diag(tanh(FFNNObj.Theta2))) * FFNNObj.W2 * (diag(tanh(FFNNObj.Theta1))); % Insert d/db1.
db2 = tanh(FFNNObj.Theta2); % Insert d/db2.
dW1 = FFNNObj.W3 * (diag(tanh(FFNNObj.Theta2))) * FFNNObj.W2 * (diag(tanh(FFNNObj.Theta1))) .* Cfg_vec; % Insert d/dW1.
dW2 = (tanh(FFNNObj.Theta2) * log(cosh(FFNNObj.Theta1)).'); % Insert d/dW2.

for a = 1:Alpha1
    dLogp(a) = sum(db1((1:Nv)+(a-1)*Nv));
end
dLogp((1:Nh2)+Alpha1) = db2;
for a = 1:Alpha1
    for b = 1:numel(BondMap)
        for n = 1:Nv
            PInd1 = Alpha1 + Nh2 + n + (a-1)*Nv;
            VInd = BondMap{b}(n); HInd = b + (a-1)*Ng;
            dLogp(PInd1) = dLogp(PInd1) + dW1(VInd,HInd);
        end
    end
    for m = 1:Nh2
        % Each 2nd layer unit has Alpha1 unique couplings.
        PInd2 = a + (m-1)*Alpha1; % Arrange by 2nd level unit index first.
        dLogp(PInd2) = sum(dW2(m,(1:Ng)+(a-1)*Ng));
    end
end

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;

% 'logderiv'