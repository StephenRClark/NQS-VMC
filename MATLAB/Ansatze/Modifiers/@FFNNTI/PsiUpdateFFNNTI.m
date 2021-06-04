% --- General FFNN wave function update function ---

function FFNNObj = PsiUpdateFFNNTI(FFNNObj,dP)
% This function updates the FFNN parameters of the ansatz from a vector of
% parameters P.
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
Alpha = FFNNObj.Alpha1;
% Unpack the changes in parameters of the NQS:
db1 = dP(1:Alpha);
db2 = dP((Alpha+1):(Alpha+Nh2));
dW1 = reshape(dP((Alpha+Nh2+1):(Alpha+Nh2+(Nv*Alpha))),Alpha,Nv);
dW2 = reshape(dP((Alpha+Nh2+(Nv*Alpha)+1):(Alpha+Nh2+(Nv*Alpha)+(Alpha*Nh2))),Nh2,Alpha);

% Apply updates to the ansatz:
FFNNObj.b1r = FFNNObj.b1 + db1;
FFNNObj.b2 = FFNNObj.b2 + db2;
FFNNObj.W1r = FFNNObj.W1 + dW1;
FFNNObj.W2r = FFNNObj.W2 + dW2;

cap = FFNNObj.ParamCap;

% Sanity check the values of the ansatz:
FFNNObj.b1r(isinf(FFNNObj.b1r)) = 0;
FFNNObj.b1r(isnan(FFNNObj.b1r)) = 0;
ind = abs(FFNNObj.b1r)>cap;
FFNNObj.b1r(ind) = sign(FFNNObj.b1r(ind))*cap;

FFNNObj.b2(isinf(FFNNObj.b2)) = 0;
FFNNObj.b2(isnan(FFNNObj.b2)) = 0;
ind = abs(FFNNObj.b2)>cap;
FFNNObj.b2(ind) = sign(FFNNObj.b2(ind))*cap;

FFNNObj.W1r(isinf(FFNNObj.W1r)) = 0;
FFNNObj.W1r(isnan(FFNNObj.W1r)) = 0;
ind = abs(FFNNObj.W1r)>cap;
FFNNObj.W1r(ind) = sign(FFNNObj.W1r(ind))*cap;

FFNNObj.W2(isinf(FFNNObj.W2)) = 0;
FFNNObj.W2(isnan(FFNNObj.W2)) = 0;
ind = abs(FFNNObj.W2)>cap;
FFNNObj.W2(ind) = sign(FFNNObj.W2(ind))*cap;

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